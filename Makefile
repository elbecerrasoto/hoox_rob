.PHONY benchmark:
benchmark: genomes.txt
	hyperfine -M 32 './hmmer100.py queries genomes.txt 100.tsv' './hmmer50.py queries genomes.txt 50.tsv'

genomes.txt:
	fd -e faa >| genomes.txt


.PHONY docs:
docs: rob_answer.pdf rob_answer.html


rob_answer.%: rob_answer.org
	pandoc rob_answer.org -o $@


.PHONY clean:
clean:
	rm rob_answer.pdf rob_answer.html genomes.txt
