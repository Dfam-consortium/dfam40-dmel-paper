# Dfam 4.0 *Drosophila melanogaster* transposable element libraries

Code for:

> Updated Transposable Element Libraries for *Drosophila melanogaster* in Dfam 4.0
>
> Clément Goubert, Anthony Gray, Robert Hubley, Travis J. Wheeler, Arian A. F. Smit

## Scripts

### `Xmatch.pl`

Runs cross_match with shorter option names and less noise, and reports the mismatch
level as mismatches over aligned bases rather than over query length. `-perbasescore`
adds a column giving the score as a percentage of the query's score against itself,
and `-cg` discounts one transition per CpG site. With `-blast` it runs rmblastn
through rmblast.pl instead of cross_match.

### `DeleteAlmostDupsfromSeqfile.pl`

Self-compares a fasta file and drops sequences that another sequence already covers,
either by divergence (`-divergence`) or by score relative to the self-score
(`-comparative`). It can also merge overlapping fragments. For DNA it calls
`Xmatch.pl` from the same directory.

Run either script with no arguments for the full option list.

## Dependencies

Perl 5, core modules only. Both scripts expect these on the PATH:

| Tool | Comes with | Needed for |
| --- | --- | --- |
| `cross_match` | phrap | everything except `-blast` |
| `rmblast.pl` | RepeatModeler, in `util/` | `-blast` |
| `rmblastn`, `makeblastdb` | rmblast | `-blast` |

`rmblast.pl` sits in RepeatModeler's `util/` subdirectory, which is often not on the
PATH:

```sh
export PATH="/usr/local/RepeatModeler/util:$PATH"
```

Each script checks for the tools its options need before writing anything, and exits
naming whatever is missing.

## Matrices

`Matrices/` holds the scoring matrices. The scripts look for it next to themselves, so
keep the two scripts and `Matrices/` in one directory.

`Matrices/crossmatch/` holds the cross_match-format matrices and `Matrices/ncbi/nt/`
the rmblastn versions of the same ones. Under `-blast`, `Xmatch.pl` swaps the
directory in a `-matrix` path, so a matrix you name on the command line has to exist
in both.

## Tested with

| | |
| --- | --- |
| cross_match | 1.090518 |
| RepeatModeler | 2.0.7 |
| rmblastn | 2.17.1+ (BLAST 2.17.0) |
| Perl | 5.34.0 |
