# Seth's GIT based version of mtsieve

Rogue maintains the official copy at https://sourceforge.net/projects/mtsieve/

To get the current version of mtsieve along with documentation for it and programs building upon it,
go to http://www.mersenneforum.org/rogue/mtsieve.html.

This repository was created in 2022-11 and it's unclear if I'll be able to continue to pull SVN changes.

It's mainly for some personal development work on core/HashTable and `sierpinski_riesel`.


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
