/*
 * SequenceSplitter.cpp
 *
 *  Created on: 20.05.2015
 *      Author: marius
 */


#include "../../include/gerbil/SequenceSplitter.h"
#include <iostream>
#include <string>
#include <fstream>


#define SS_MMER_UNDEFINED 0x80000000
#define SS_BIN_UNDEFINED 0xffff
#define SS_MMER_VAL_UNDEFINED 0x80000000

gerbil::SequenceSplitter::SequenceSplitter(
		const uint32 &superBundlesNumber,					// number of SUperBundles
		SyncSwapQueueMPMC<ReadBundle>* readBundleSyncQueue,	// Queue with Bundles of reads
		const uint8 &splitterThreadsNumber,					// Number of SplitterThreads
		const uint32_t &k,									// size of k-mer
		const uint8 &m,										// Size of m-mer (minimizer)
		const uint_tfn &tempFilesNumber,					// Number of tempFiles
		const bool &norm
): 	_norm(norm), _k(k), _m(m), _tempFilesNumber(tempFilesNumber), _splitterThreadCount(splitterThreadsNumber),
	_readBundleSyncQueue(readBundleSyncQueue), _leftSuperBundles(NULL), _baseNumbers(0)
{
	_superBundleQueues = new SyncSwapQueueMPSC<SuperBundle>*[SB_WRITER_THREADS_NUMBER];
		for(uint i = 0; i < SB_WRITER_THREADS_NUMBER; ++i)
			_superBundleQueues[i] = new SyncSwapQueueMPSC<SuperBundle>(superBundlesNumber / SB_WRITER_THREADS_NUMBER);

	_splitterThreads = new std::thread*[_splitterThreadCount];

	uint32 mMerNumber = 1 << (2 *_m); 	// 4^m

	_mVal = new uint32[mMerNumber];
	_mToBin = new uint_tfn[mMerNumber];
}

gerbil::SequenceSplitter::~SequenceSplitter() {
	for(uint i = 0; i < SB_WRITER_THREADS_NUMBER; ++i)
		delete _superBundleQueues[i];
	delete[] _superBundleQueues;
	delete[] _splitterThreads;
	delete[] _mVal;
	delete[] _mToBin;
}


gerbil::SyncSwapQueueMPSC<gerbil::SuperBundle>** gerbil::SequenceSplitter::getSuperBundleQueues() {
	return _superBundleQueues;
}
uint32_t gerbil::SequenceSplitter::invMMer(const uint32_t &mmer){
	uint32 rev = 0;
	uint32 immer = ~mmer;
	for(uint32 i = 0 ; i < _m ; ++i)
	{
		rev <<= 2;
		rev |= immer & 0x3;
		immer >>= 2;
	}
	return rev;
}


// kmc2 strategy
bool gerbil::SequenceSplitter::isAllowed(uint32 mmer) {
	if(mmer > invMMer(mmer))
		return false;

	if(mmer == (0xaaaaaaaa >> (2 * (16 - _m + 2)))) 	// AAGGG..G
		return false;
	if(mmer== (0x55555555 >> (2 * (16 - _m + 2)))) 		// AACCC..C
		return false;

	for (uint32 j = 0; j < _m - 3; ++j)
		if ((mmer & 0xf) == 0) 							// AA inside
			return false;
		else
			mmer >>= 2;

	if ((mmer & 0xf) == 0)    							// *AA prefix
		return false;
	if (mmer == 0)           			 				// AAA prefix
		return false;
	if (mmer == 0x04)        							// ACA prefix
		return false;

	return true;
}

void gerbil::SequenceSplitter::loadUhs(std::string path)
{
	uint32 mmersNumber = 1<<(2*_m);
	uhs = new byte[mmersNumber];
	for(int i = 0; i< mmersNumber; i++) 
    	uhs[i] = 0;
    
    char num_codes[256];
	for (int i = 0; i < 256; i++)
        num_codes[i] = -1;
	num_codes['A'] = num_codes['a'] = 0;
	num_codes['C'] = num_codes['c'] = 1;
	num_codes['G'] = num_codes['g'] = 2;
	num_codes['T'] = num_codes['t'] = 3;

	std::ifstream infile(path);
	int str_value;
	std::string line;
	std::cout << "lolz0" << std::endl<< std::flush;
	while (std::getline(infile, line))
	{
		str_value = 0;
		const char *seq = line.c_str();
		for(uint32 i = 0; i < _m; i++)
		{
			str_value += num_codes[seq[i]] << (2*_m - (1+i)*2);
		}
		//str_value = (num_codes[seq[0]] << 16) + (num_codes[seq[1]] << 14) + (num_codes[seq[2]] << 12) + (num_codes[seq[3]] << 10) + (num_codes[seq[4]] << 8) + (num_codes[seq[5]] << 6) + (num_codes[seq[6]] << 4) + (num_codes[seq[7]]<< 2) + (num_codes[seq[8]]);
		uhs[str_value] = 1;
		uhs[invMMer(str_value)] = 1;
	}
	std::cout << "lolz6" << std::endl<< std::flush;
	infile.close();
	std::cout << "lolz7" << std::endl<< std::flush;
}

void gerbil::SequenceSplitter::loadRanks(uint32* ranks)
{
	uint32 mmersNumber = 1<<(2*_m);
	uint32 localRanks[mmersNumber];
	std::ifstream infile("ranks.txt");
	int str_value;
	std::string line;
	uint32 i = 0;
	while (std::getline(infile, line))
	{
		str_value = std::stoi(line);
		ranks[i] = str_value;
		localRanks[i] = i;
		i++;
		if(i == mmersNumber) break;
	}
	std::sort(localRanks, localRanks + mmersNumber, CompR(ranks));
	for(uint32 j = 0; j < mmersNumber; ++j)
	{
		std::cout << " j is = " << j << " val is " << localRanks[j] << std::endl;
		ranks[localRanks[j]] = j;
	}
	infile.close();
}

void gerbil::SequenceSplitter::loadFrequency(uint32* freq)
{
	uint32 mmersNumber = 1<<(2*_m);
	std::ifstream infile("freq.txt");
	int str_value;
	std::string line;
	uint32 i = 0;
	while (std::getline(infile, line))
	{
		str_value = std::stoi(line);
		freq[i] = str_value;
		i++;
		if(i == mmersNumber) break;
	}
	infile.close();
}


void gerbil::SequenceSplitter::detMMerHisto() {
	const uint32 mMersNumber = 1 << (2 * _m);
	uint32 stats[mMersNumber];

	loadUhs();
	loadRanks(_mVal);
	
	loadFrequency(stats);

	// kmc2 strategy
	// std::vector<std::pair<uint32_t,uint32_t>> kMerFrequencies;
	// for(uint32 mmer = 0; mmer < mMersNumber; ++mmer)
	// 	kMerFrequencies.push_back(std::make_pair(isAllowed(mmer) ? mmer : 0xffffffffu, mmer));
	//  uint32 normalized_value;
	//  for(uint32 mmer = 0; mmer < mMersNumber; ++mmer)
	//  {
	//  	if(uhs[mmer] == 0)
	//  	{
	//  		if(isAllowed(mmer))
	//  			normalized_value = 2*mMersNumber + mmer;
	//  		else
	//  			normalized_value = 3*mMersNumber + mmer;
	//  	}
	//  	else
	//  	{
	//  		if(isAllowed(mmer))
	//  			normalized_value = mmer;
	//  		else
	//  			normalized_value = mMersNumber + mmer;
	//  	}
	//  	kMerFrequencies.push_back(std::make_pair(normalized_value, mmer));	
	//  }
		

	// std::sort(kMerFrequencies.begin(), kMerFrequencies.end());
	// uint32 rank = 0;
	// for (std::vector<std::pair<uint32_t,uint32_t>>::iterator it = kMerFrequencies.begin() ; it != kMerFrequencies.end(); ++it) {
	// 	_mVal[it->second] = rank++;
	// }

	

	// for(uint32 i = 0; i < mMersNumber; ++i)
	// {
	// 	_mToBin[i] = i % _tempFilesNumber;
	// }


	// pyramid-shaped distribution
	// uint32 bfnx2 = 2 * _tempFilesNumber;
	// uint32 mMN = mMersNumber - mMersNumber % bfnx2;
	// int32 p = mMN / bfnx2;
	// for(uint32 i = 0; i < mMersNumber; ++i)
	// 	if(_mVal[i] < mMN){
	// 		int32 d = _mVal[i] / bfnx2;
	// 		int32 m = (_mVal[i] + d * bfnx2 / p) % bfnx2;
	// 		int32 x = m - (int32)_tempFilesNumber;
	// 		// if((x < 0 ? -x - 1 : x) == 335) std::cout << i << std::endl;
	// 		_mToBin[i] = x < 0 ? -x - 1 : x;
	// 	}
	// 	else
	// 	{
	// 		// if((_tempFilesNumber - 1 - (_mVal[i] % _tempFilesNumber)) == 335) std::cout << i << std::endl;
	// 		_mToBin[i] = _tempFilesNumber - 1 - (_mVal[i] % _tempFilesNumber);
	// 	}


	uint32 local_sorted[mMersNumber];
	uint32 *sorted = local_sorted;
	for (uint32 i = 0; i < mMersNumber ; ++i)
		sorted[i] = i;
	std::sort(sorted, sorted + mMersNumber, Comp(stats));

	std::list<std::pair<uint32, uint64>> _stats;

	for (uint32 i = 0; i < mMersNumber ; ++i)
	{
		_stats.push_back(std::make_pair(sorted[i], stats[sorted[i]]));
	}

	std::list<uint32> _stats_zero_ids;
	for (auto it = _stats.begin(); it != _stats.end();)
	{
		if (it->second == 0)
		{
			_stats_zero_ids.push_back(it->first);
			it = _stats.erase(it);
		}
		else
			++it;
	}
	uint32 total_zeroes = _stats_zero_ids.size();
	uint32 assigned_zeros = 0;
	

	std::list<std::pair<uint32, uint64>> group;
	uint32 bin_no = 0;
	//counting sum
	double sum = 0.0;
	for (auto &i : _stats)
	{
		// i.second += 1000;
		//i.second += 100;
		sum += i.second;
	}


	double mean = sum / _tempFilesNumber;
	double max_bin_size = mean;
	uint32 n = _tempFilesNumber; //one is needed for disabled signatures
	uint32 max_bins = _tempFilesNumber;
	uint32 num_zeros_in_bin = (_stats_zero_ids.size() / max_bins) + 1;
	
	//max_bin_size = (mean + _stats.front().second) / 2;

	while (_stats.size() > n)
	{
		std::pair<uint32, uint64>& max = _stats.front();

		if (max.second >= mean)
		{
			std::cout << max.first << " innn " << bin_no <<std::endl << std::flush;
			_mToBin[max.first] = bin_no++;
			sum -= max.second;
			mean = sum / (max_bins - bin_no);
			//max_bin_size = 1.1 * mean;
			//max_bin_size = _stats.back().second;

			_stats.pop_front();
			--n;
		}
		else	
		{
			//heuristic
			group.clear();
			double tmp_sum = 0.0;
			uint32 in_current = 0;
			for (auto it = _stats.begin(); it != _stats.end();)
			{
				if (tmp_sum + it->second <= max_bin_size)
				{
					tmp_sum += it->second;
					group.push_back(*it);
					it = _stats.erase(it);
					++in_current;
				}
				else
					++it;
			}
			for (auto i = group.begin(); i != group.end(); ++i)
			{
				//std::cout << i->first << " in " << bin_no <<std::endl << std::flush;
				_mToBin[i->first] = bin_no;
			}

			for (auto p = std::make_pair(_stats_zero_ids.begin(),0);
			 	p.first != _stats_zero_ids.end() && p.second < num_zeros_in_bin && assigned_zeros < total_zeroes;
			  	++p.second, ++assigned_zeros)
			{
				_mToBin[*(p.first)] = bin_no;
				p.first = _stats_zero_ids.erase(p.first);
			}
			--n;
			++bin_no;

			sum -= tmp_sum;
			mean = sum / (max_bins - bin_no);
			//if(bin_no > _tempFilesNumber / 4)
				max_bin_size = 1.01 * mean;
			//if(bin_no < 200)
			//max_bin_size = _stats.back().second;
			
		}
	}
	if (_stats.size() > 0)
	{
		for (auto i = _stats.begin(); i != _stats.end(); ++i)
		{
			_mToBin[i->first] = _tempFilesNumber - 1 - (i->first%_tempFilesNumber);
		}
	}
	if(_stats_zero_ids.size() == 0)
	{
		std::cout << "No elements in zero sized " << std::endl;
	}
}

void gerbil::SequenceSplitter::process() {
	// calculate statistics
	detMMerHisto();

	_leftSuperBundles = new SuperBundle**[_tempFilesNumber];
	for(uint i(0); i < _tempFilesNumber; ++i)
		_leftSuperBundles[i] = new SuperBundle*[_splitterThreadCount];

	// split reads to s-mers
	for(uint i = 0; i < _splitterThreadCount; ++i)
		_splitterThreads[i] = new std::thread([this](uint id){
			processThread(id);
		}, i);
}

void gerbil::SequenceSplitter::join() {
	for(uint i = 0; i < _splitterThreadCount; ++i) {
		_splitterThreads[i]->join();
		delete _splitterThreads[i];
	}

	SuperBundle* lastSuperBundle;

	// merge left superbundles
	for(uint_tfn tempFileId(0); tempFileId < _tempFilesNumber; ++tempFileId) {
		lastSuperBundle = _leftSuperBundles[tempFileId][0];
		// merge all superbundles for unique file
		for(size_t tId(1); tId < _splitterThreadCount; ++tId) {
			// copy low to big
			if(lastSuperBundle->getSize() < _leftSuperBundles[tempFileId][tId]->getSize())
				std::swap(lastSuperBundle, _leftSuperBundles[tempFileId][tId]);
			// merge
			if(!lastSuperBundle->merge(*_leftSuperBundles[tempFileId][tId])) {
				// store at fail
				lastSuperBundle->finalize();
				if(!lastSuperBundle->isEmpty())
					_superBundleQueues[tempFileId % SB_WRITER_THREADS_NUMBER]->swapPush(lastSuperBundle);
				std::swap(lastSuperBundle, _leftSuperBundles[tempFileId][tId]);
			}
		}
		// store last superbundle
		lastSuperBundle->finalize();
		if(!lastSuperBundle->isEmpty())
			_superBundleQueues[tempFileId % SB_WRITER_THREADS_NUMBER]->swapPush(lastSuperBundle);
	}

	// finalize all queues
	for(uint i = 0; i < SB_WRITER_THREADS_NUMBER; ++i)
		_superBundleQueues[i]->finalize();

	// free memory
	for(uint_tfn tempFileId(0); tempFileId < _tempFilesNumber; ++tempFileId) {
		for(size_t tId(1); tId < _splitterThreadCount; ++tId)
			delete _leftSuperBundles[tempFileId][tId];
		delete[] _leftSuperBundles[tempFileId];
	}
	delete[] _leftSuperBundles;
}

void gerbil::SequenceSplitter::print() {
	printf("number of bases        : %12lu\n", _baseNumbers.load());
}

//private

/*
 * adds a s-mer to the corresponding SuperBundle
 * if necessary, the SuperBundle is stored in the queue
 */
#define SS_SAVE_S_MER {																								\
	curTempFileId = bins[min_val_pos];																				\
	if(!curSuperBundles[curTempFileId]->add(rb->data + i - smer_c - _k + 1, smer_c + _k - 1, _k)) {					\
		curSuperBundles[curTempFileId]->finalize();																	\
		IF_MESS_SEQUENCESPLITTER(sw.hold();)																		\
		_superBundleQueues[curTempFileId % SB_WRITER_THREADS_NUMBER]->swapPush(curSuperBundles[curTempFileId]);		\
		IF_MESS_SEQUENCESPLITTER(sw.proceed();)																		\
		curSuperBundles[curTempFileId]->tempFileId = curTempFileId;													\
		curSuperBundles[curTempFileId]->add(rb->data + i - smer_c - _k + 1, smer_c + _k - 1, _k);					\
	}}

/*
 * 0x0 => A, a
 * 0x1 => C, c
 * 0x2 => G, g
 * 0x3 => T, t
 * 0x4 => 'N', 'n', '.'
 * 0x5 => 'E' (end of line)
 * 0xf => others (Error)
 */
const gerbil::byte baseToByteA[256] = {
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0x4, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0x0, 0xf, 0x1, 0xf, 0x5, 0xf, 0x2, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0x4, 0xf,
		0xf, 0xf, 0xf, 0xf, 0x3, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0x0, 0xf, 0x1, 0xf, 0xf, 0xf, 0x2, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0x4, 0xf,
		0xf, 0xf, 0xf, 0xf, 0x3, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
		0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf
};

void gerbil::SequenceSplitter::processThread(const uint &id) {
	IF_MESS_SEQUENCESPLITTER(
		StopWatch sw(CLOCK_THREAD_CPUTIME_ID);
		sw.start();
	)
	ReadBundle* rb = new ReadBundle;
	uint32 rc;

	byte* reads_p_end;
	byte* cur_p;
	byte cur;
	uint32 mmer;
	uint32 inv_mmer;

	uint32 m_mask = 0;
	for(uint i = 0; i < _m; ++i)
		m_mask = (m_mask << 2) | 0x3;

	const byte inv[4] = {0x3, 0x2, 0x1, 0x0};

	uint32* bins = new uint32[READ_BUNDLE_SIZE_B];
	uint32* bins_p;
	uint32* mmerval = new uint32[READ_BUNDLE_SIZE_B];
	uint32* mmerval_p;
	uint32 reads_size;
	uint32 mmer_min;

	uint64 baseNumbers(0);


	byte ic;

	SuperBundle** curSuperBundles = new SuperBundle*[_tempFilesNumber];
	for(uint_tfn i = 0; i < _tempFilesNumber; ++i) {
		curSuperBundles[i] = new SuperBundle();
		curSuperBundles[i]->tempFileId = i;
	}

	while(_readBundleSyncQueue->swapPop(rb)) {
#ifdef DEB_SS_LOG_RB
		rb->print();
#endif
		rc = *(rb->readsCount);
		reads_size = *(rb->readOffsets - rc);
		reads_p_end = rb->data + reads_size;
		cur_p = rb->data;
		ic = _m;
		mmer = 0;
		inv_mmer = 0;
		mmerval_p = mmerval;
		bins_p = bins;
		// save mmer values
		for(; cur_p < reads_p_end; ++cur_p, ++mmerval_p, ++bins_p) {
			*cur_p = (cur = baseToByteA[*cur_p]);
			if(ic)
				--ic;
			if(cur < (byte)0x4) {
				mmer = ((mmer << 2) | cur) & m_mask;
				inv_mmer = (inv_mmer >> 2) | (inv[cur] << (2*(_m-1)));
				++baseNumbers;
			}
			else {
				ic = _m;
			}
			if(ic) {
				*mmerval_p = SS_MMER_VAL_UNDEFINED;
				*bins_p = SS_BIN_UNDEFINED;
			}
			else {
				mmer_min = (!_norm || mmer < inv_mmer) ? mmer : inv_mmer;
				*mmerval_p = _mVal[mmer_min];
				*bins_p = _mToBin[mmer_min];
			}

		}

		//
		uint32 min_val = SS_MMER_VAL_UNDEFINED;
		uint32 min_val_pos = 0;
		uint cur_val;
		int32 smer_c = (uint32) _m - _k;
		uint_tfn curTempFileId;

		for(uint32 i(_m - 1); i < reads_size; ++i) {
			cur_val = mmerval[i];


			if(cur_val == SS_MMER_VAL_UNDEFINED) {
				// completed s-mer --> undefined minimizer (e.g. end of read)
				// save the last smer
				if(smer_c > 0) {
					SS_SAVE_S_MER;
				}
				smer_c = _m - _k;
				if(rb->data[i] > (byte)0x3) {
					// jump
					i += _m - 1;
				}

				min_val = SS_MMER_VAL_UNDEFINED;
			}
			else if(min_val == SS_MMER_VAL_UNDEFINED) {
				// beginning of a new s-mer
				min_val = cur_val;
				min_val_pos = i;
				++smer_c;
			}
			else if(min_val == cur_val) {
				// extended s-mer --> equal minimizer
				min_val_pos = i;
				++smer_c;
			}
			else if(cur_val < min_val) {
				// extended s-mer --> better minimizer, but same bin
				// completed s-mer --> better minimizer
				if(smer_c > 0 && bins[min_val_pos] != bins[i]) {
					// save the last s-mer
					SS_SAVE_S_MER;
					smer_c = 1;
				}
				else
					++smer_c;
				min_val = cur_val;
				min_val_pos = i;
			}
			else if(i - min_val_pos > _k - _m || smer_c > 512) {
			//else if(i - min_val_pos > _k - _m) {
				// completed s-mer --> minimizer is out of range
				// save last s-mer
				if(smer_c > 0) {
					SS_SAVE_S_MER;
					smer_c = 1;
				}
				else
					++smer_c;
				min_val = cur_val;
				min_val_pos = i;
				for(uint x = _k - _m; x > 0 ; --x)
					if(mmerval[i - x] < min_val) {
						min_val = mmerval[i - x];
						min_val_pos = i - x;
					}
			}
			else {
				// extended s-mer --> current minimizer is still the best
				++smer_c;
			}
		}
		rb->clear();
	}

	_baseNumbers += baseNumbers;

	for(uint_tfn i = 0; i < _tempFilesNumber; ++i)
		_leftSuperBundles[i][id] = curSuperBundles[i];

	// free memory
	delete rb;
	delete[] bins;
	delete[] mmerval;
	delete[] curSuperBundles;


	IF_MESS_SEQUENCESPLITTER(
		sw.stop();
		printf("splitter[%2u]: %.3f s\n", id, sw.get_s());
	)

}
