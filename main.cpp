#include "kmer.h"
#include "sequence_container.h"
#include <iostream>
#include <unordered_set>
#include <map>


int main(int argc, char** argv)
{
	const int32_t WINDOW = 100;
	const int32_t MIN_PART = 500;
	const int32_t MAX_CUTS = 2;

	if (argc < 2)
	{
		std::cerr << "Usage: pbclip reads_1 [reads_2 ...]\n";
		return 1;
	}

	std::vector<std::string> readsList;
	for (int i = 1; i < argc; ++i)
	{
		readsList.emplace_back(argv[i]);
	}

	Parameters::get().numThreads = 1;
	Parameters::get().kmerSize = 15;
	Parameters::get().minimumOverlap = 0;
	Parameters::get().unevenCoverage = false;

	SequenceContainer readsContainer;
	try
	{
		for (auto& readsFile : readsList)
		{
			readsContainer.loadFromFile(readsFile, /*min length*/ 0);
		}
	}
	catch (SequenceContainer::ParseException& e)
	{
		std::cerr << e.what();
		return 1;
	}
	std::cerr << "Sequences loaded\n";

	std::unordered_map<Kmer, size_t> kmerConunter;
	//std::unordered_set<Kmer> usedKmers;
	std::vector<int32_t> splitPos;
	std::map<int32_t, int32_t> hist;
	int numGood = 0;
	int numChopped = 0;
	int numComplex = 0;
	std::vector<FastaRecord> outputSequences;
	for (auto& seq : readsContainer.iterSeqs())
	{
		if (!seq.id.strand()) continue;

		kmerConunter.clear();
		splitPos.clear();
		//usedKmers.clear();
		hist.clear();
		size_t duplicates = 0;
		for (const auto& kmerPos : IterKmers(seq.sequence))
		{
			Kmer stdKmer(kmerPos.kmer);
			stdKmer.standardForm();

			if (kmerConunter.count(kmerPos.kmer))
			{
				++duplicates;
				//if (!usedKmers.count(stdKmer))
				{
					splitPos.push_back((kmerPos.position + 
									    kmerConunter[kmerPos.kmer]) / 2);
				}
				//usedKmers.insert(stdKmer);
			}

			kmerConunter[kmerPos.kmer.reverseComplement()] = kmerPos.position;
		}

		for (int32_t x: splitPos) ++hist[x / WINDOW];

		//std::cerr << seq.description << " " << seq.sequence.length() 
		//	<< " " << duplicates << std::endl;

		/*for (auto x: splitPos)
		{
			std::cerr << x << " ";
		}
		std::cerr << std::endl;*/

		std::vector<int32_t> chopPoints;
		for (auto window : hist)
		{
			if (window.second > 25)
			{
				int32_t pos = window.first * WINDOW;
				int32_t prevChop = chopPoints.empty() ? 0 : chopPoints.back();
				if (pos > prevChop + MIN_PART && 
					pos < (int)seq.sequence.length() - MIN_PART)
				{
					chopPoints.push_back(pos);
					//std::cerr << window.second << std::endl;
				}

				//std::cerr << "Anomaly: " << window.first * WINDOW 
				//	<< " " << window.second << std::endl;
			}
		}

		/*for (auto x: chopPoints)
		{
			std::cerr << "Chop: "<< x << std::endl;
		}*/

		if (chopPoints.empty()) 
		{
			++numGood;
			outputSequences.push_back(seq);
		}
		else if (chopPoints.size() < MAX_CUTS) 
		{
			++numChopped;
			chopPoints.push_back(seq.sequence.length());
			for (size_t i = 0; i < chopPoints.size(); ++i)
			{
				int32_t start = i > 0 ? chopPoints[i - 1] : 0;
				DnaSequence subread = seq.sequence.substr(start, chopPoints[i] - start);
				FastaRecord rec(subread, seq.description + 
								"_subread_" + std::to_string(i), FastaRecord::ID_NONE);
				outputSequences.push_back(rec);
			}
		}
		else 
		{
			++numComplex;
		}
	}

	std::cerr << "Good: " << numGood << " chopped: " << numChopped
		<< " bad: " << numComplex << std::endl;
	SequenceContainer::writeFasta(outputSequences, "");

	return 0;
}
