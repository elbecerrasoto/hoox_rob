.PHONY benchmark:
benchmark: genomes.txt
	hyperfine -M 16 './hmmer100.py queries genomes.txt 100.tsv' './hmmer50.py queries genomes.txt 50.tsv'


.PHONY benchmark-few:
benchmark-few: genomes.txt
	hyperfine -M 16 './hmmer100.py queries_few genomes.txt 100.tsv' './hmmer50.py queries_few genomes.txt 50.tsv'


genomes.txt:
	fd -I -e faa >| genomes.txt


.PHONY docs:
docs: rob_answer.pdf rob_answer.html


rob_answer.%: rob_answer.org
	pandoc rob_answer.org -o $@


.PHONY clean:
clean:
	rm -f rob_answer.pdf rob_answer.html genomes.txt 50.tsv 100.tsv
