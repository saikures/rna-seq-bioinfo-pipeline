import sys

# Usage: python prepare_links.py input.paf > links.txt

with open(sys.argv[1]) as paf:
    for line in paf:
        parts = line.strip().split('\t')
        query, qlen, qstart, qend, strand, ref, rlen, rstart, rend = parts[:9]
        print(f"{ref}\t{rstart}\t{rend}\t{query}\t{qstart}\t{qend}")
