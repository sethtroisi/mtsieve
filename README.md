# Seth's GIT based version of mtsieve

Rogue maintains the official copy at https://sourceforge.net/projects/mtsieve/

To get the current version of mtsieve along with documentation for it and programs building upon it,
go to http://www.mersenneforum.org/rogue/mtsieve.html.

This repository was created in 2022-11 and it's unclear if I'll be able to continue to pull SVN changes.

It's mainly for some personal development work on core/HashTable and `sierpinski_riesel`.


### Testing

```
# $ cat crus_sequences.txt  | grep '^[0-9]' | shuf -n 1 -o crus_seqs_rand1.txt
# $ cat crus_sequences.txt  | grep '^[0-9]' | shuf -n 10 -o crus_seqs_rand10.txt
# $ cat crus_sequences.txt  | grep '^[0-9]' | shuf -n 100 -o crus_seqs_rand100.txt
# $ cat crus_sequences.txt  | grep '^[0-9]' | shuf -n 1000 -o crus_seqs_rand1000.txt

$ wc crus_seqs_rand*
 1000  3000 15978 crus_seqs_rand1000.txt
  100   300  1534 crus_seqs_rand100.txt
   10    30   158 crus_seqs_rand10.txt
    1     3    13 crus_seqs_rand1.txt

rm temp*.out
for N in {1,10,100,1000}; do
  cat "crus_seqs_rand${N}.txt" | awk -F", " '{ print "-s\"" $1 "*" $2 "^n" $3 "\"" }' | tr '\n' ' ' > "seqs${N}.txt"
  wc "seqs${N}.txt"
  echo
  eval "./srsieve2_clean -P 1e8 -W8 -n1e5 -N2e5 $(cat seqs${N}.txt) -o temp_${N}.in"
  echo
  eval "./srsieve2_clean -P 2e8     -i temp_${N}.in -o temp_${N}_befor1.out"
  echo
  eval "./srsieve2_clean -P 1e9 -W8 -i temp_${N}.in -o temp_${N}_befor8.out"
  echo
  eval "./srsieve2       -P 2e8     -i temp_${N}.in -o temp_${N}_after1.out"
  echo
  eval "./srsieve2       -P 1e9 -W8 -i temp_${N}.in -o temp_${N}_after8.out"
  echo -e "\n\n"
done

```

### primesieve

Kim Walisch's amazing primesieve is included with a submodule.

You may need to `git submodule update --init` after cloning this repository

To match the old structure I had to do this (you don't have to)

```
mv sieve sieve_old

# Old approach that leaves primesieve dirty
ln -s submodules/primesieve/src/ sieve
ln -s ../include/primesieve.{h,hpp} submodules/primesieve/src/

# New approach with more symlinks
mkdir sieve
cd sieve
ln -s ../submodules/primesieve/src/*.cpp .
ln -s ../submodules/primesieve/include/* .
ln -s ../submodules/primesieve/include/primesieve .
```
