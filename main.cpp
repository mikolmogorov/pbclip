//(c) 2016 by Authors
//Released under the BSD license (see LICENSE file)

#include "kmer.h"
#include "sequence_container.h"
#include <iostream>
#include <unordered_set>
#include <map>




int main(int argc, char** argv)
{
	//const int32_t WINDOW = 100;
	const int32_t MIN_PART = 500;
	const int32_t MAX_CUTS = 2;
	const int MIN_SUSPICIOUS_KMERS = 10;
	const int LINKAGE_DIST = 50;
	const int MAX_CHOP_RANGE = 100;

	if (argc < 2)
	{
		std::cerr << "Usage: pbclip reads_file_1 [reads_file_2 ...]\n\n"
			<< "Fix PacBio reads that were not properly split into subreads.\n"
			<< "Input reads should be in fasta/q format, could be gzipped.\n"
			<< "Outputs fasta to stdout.";
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
		hist.clear();
		for (const auto& kmerPos : IterKmers(seq.sequence))
		{
			Kmer stdKmer(kmerPos.kmer);
			stdKmer.standardForm();

			if (kmerConunter.count(kmerPos.kmer))
			{
				splitPos.push_back((kmerPos.position + 
									kmerConunter[kmerPos.kmer]) / 2);
			}

			kmerConunter[kmerPos.kmer.reverseComplement()] = kmerPos.position;
		}

		std::sort(splitPos.begin(), splitPos.end());

		/*std::cerr << seq.description << " " << seq.sequence.length() 
			<< " " << splitPos.size() << std::endl;
		for (auto x: splitPos)
		{
			std::cerr << x << " ";
		}
		std::cerr << std::endl;*/

		splitPos.push_back(seq.sequence.length() - 1);
		std::vector<int32_t> chopPoints;
		int lastEnd = 0;
		int lastStart = 0;
		int clustSize = 0;
		int prevChop = 0;

		for (int32_t pos : splitPos)
		{
			if (pos - lastEnd < LINKAGE_DIST)
			{
				++clustSize;
				lastEnd = pos;
			}
			else
			{
				int32_t chopPos = (lastEnd + lastStart) / 2;
				if (clustSize >= MIN_SUSPICIOUS_KMERS &&
					chopPos - prevChop > MIN_PART &&
					lastEnd - lastStart < MAX_CHOP_RANGE)
				{
					chopPoints.push_back(chopPos);
					prevChop = chopPos;
					//std::cerr << "Chop: " << chopPos << std::endl;
				}
				lastStart = pos;
				lastEnd = pos;
				clustSize = 1;
			}
		}

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
	SequenceContainer::writeFasta(outputSequences, "", 
								  /*only pos strand*/ true);

	return 0;
}
