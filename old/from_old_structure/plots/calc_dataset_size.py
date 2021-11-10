from Bio import SeqIO

#Get dataset sizes in base pairs for calculating relative abundance score

if __name__ == "__main__":
    fastqs = ["../deduplication/deduplicated_sequel-demultiplex.1896_A01.ccs.fastq"
             ,"../deduplication/deduplicated_sequel-demultiplex.1896_B01.ccs.fastq"
             ,"../deduplication/deduplicated_sequel-demultiplex.1896_C01.ccs.fastq"
             ,"../deduplication/deduplicated_sequel-demultiplex.MOCK_D01.ccs.fastq"
             ,"../deduplication/deduplicated_sequel-demultiplex.MOCK_E01.ccs.fastq"
             ,"../deduplication/deduplicated_sequel-demultiplex.MOCK_F01.ccs.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152366_1.part-01.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152366_1.part-02.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152366_1.part-03.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152366_1.part-04.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-01.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-02.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-03.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-04.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-05.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-06.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-07.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-08.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-09.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-10.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-11.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-12.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-13.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-14.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-15.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-16.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-17.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-18.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-19.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-20.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-21.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-22.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-23.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-24.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-25.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-26.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-27.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-28.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-29.fastq"
             ,"/home/noyes046/jsettle/loman_mock/ERR3152367_1.part-30.fastq"]

    for fastq in fastqs:
        bp = 0
        for rec in SeqIO.parse(fastq, "fastq"):
            bp += len(rec.seq)

        print(fastq + " size (bp): " + str(float(bp)/1000000000.0))
