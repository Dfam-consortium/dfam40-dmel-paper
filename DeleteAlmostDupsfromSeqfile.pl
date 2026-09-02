#!/usr/bin/perl -w
use Getopt::Long;
use strict;
use FindBin;

# Xmatch.pl and the Matrices/ directory sit next to this script.
# cross_match, rmblast.pl, rmblastn and makeblastdb have to be on the PATH.
my $xmatch = "$FindBin::RealBin/Xmatch.pl";
my $matrixdir = "$FindBin::RealBin/Matrices/crossmatch";

#200924 added -j and -table option
#201111 fixed a bug causing multiple overlapping fragments to produce only a combination of the last two fragments
#201214 name of exact duplicate of removed sequence now in output
#210107 script now changes from the for-DNA-default Nscores1matrix to matrix when N strings of 50 or longer are in the query.
#210726 Fragments of different orientation are now merged (frags in reverse orientation to the contig are flipped first)
#210414 Need to handle more fragment descriptions
## One with :(\d+)-(\d+) where reverse strand is indicated by  $2 being smaller than $1
## Also :(\d+)-(\d+):(\d+)-(\d+) (with above rule

#211031: Fixed the line $seq =~ s/^>(\S+)\n// to $seq =~ #s/^>(\S+).+\n// at the place where exact same sequenes are compared. Ouch.
#        Kept information what is removed in non-cross_match round represented in the -r (reanalysis) standard out
#211119: Added -comparative mode (requiring access to my Xmatch.pl). Sequences are removed if they score > set % of their score-to-self
#        aligned to another sequence (restrictions of containment in other sequence are the same).
#211130: Added -keep option
#211203: Using symmetricNscores1matrix by default instead of Nscores1matrix, which assumes a42% GC-richess and a preference of GC->AT subs
#211203: Added relpairsc hash (as relative scores are not the same for A vs B and B vs A);+ stop trying once below mimimum relscore
#220505: Keeplist only worked if lines started with > or with fasta file; fixed.
#220505: replaced individual lines with removeone subroutine
#220727: perfect matches in reverse were not processed; fixed
#220805: filtering information output to <file>.filterinfo
#220810: Now sorting by score/aligned-length-of-query when not -comparative
#220820+: Now verifies that each sequence A removing sequence B also would have removed those sequences previously removed by B (otherwise B is kept)
#230915: fixed /u1/home -> /hom directories and >& redirects for new system
#231019: allowed reanalysis with a higher cutoff score
#231019: fixed horrid merging failure ( was ($seq{$secname} = $entry{$secname}) =~ s/^>.*\n([A-Za-z]+\n)/$1/;
#        now: ($seq{$secname} = $entry{$secname}) =~ s/^>.*\n([A-Za-z]+)/$1/; I probably changed the nature of %seq once
#260902: calls Xmatch.pl and reads matrices next to this script instead of
#        /home/asmit; cross_match and rmblast.pl now have to be on the PATH;
#        dropped commented-out code
#241004: included maximum bp left on both sides combined option
#        also a few fixes (e.g. at someplaces forgotten that begin query is 1-based)

## Code to check for local regions of high divergence should be included Check
## /home/asmit/Davidslibrary/Newlib/Cleanup/newchoosedoubles.pl lines
## 177-246. Should be applied with or without -comparative This involves a
## realignment of fragments of sequence 1 to full seq2 If alignments are already
## available, it would be much faster to simply calculated divergence and
## relative score from that, but that also would require writing more code...

## to do
## = option to keep mixed upper and lower case in output (requires a copy of %entry)
## - when entries are reverse identical or happy to be deleted, check which one
##   matches the most/best in forward direction to other sequences
## - in comparative mode, it may make more sense always to show in the info file
##   the cross_match line where the first sequence (query) gets dumped for the
##   secon (target), since one than sees the relative score. It also makes the
##   lines easier to interpret, even when not in comparative mode.
## Currently merging of fragments fails when there is a mix of conventions
## (sequences labeled chr1_1_200_C with chr1:101_300_R)
## There is also an error report 
## Use of uninitialized value in numeric comparison (<=>) at /home/asmit/bin/DeleteAlmostDupsfromSeqfile.pl line 622, <IN> line 294.
## Even when things go right

my $USAGE = "Usage: DeleteAlmostDupsfromSeqfile <fasta seqfile> 

Output: 
<f>.selection   Fasta file with remaining sequences
<f>.DADinfo     File listing changes to input file resulting in <f>.selection
<f>.DADmod      A modified fasta file with fixed same names, perfect duplicate
                seqs removed and overlapping sequences merged (if option -me)
<f>.DADmod.mods File listing modifications to input file resulting in <f>.DADmod
<f>.vse         A cross_match formatted self comparison of <file> or (if different) 

Options
 -p(rotein)       (default DNA)

 -d(ivergence) #  Maximum divergence in percent for elimination (default 10%)

 -b(last)         Run with rmblastn instead of cross_match. Much faster, but not yet
                  working with the -comparative option
 -c(omparative) # Minimum fraction (in %) of score to itself for elimination.
		  Alternative to -d only; only tested on DNA
                  With this option non-complexity adjusted scores are returned,
                  so the -score option may have to be set higher 
                  By default \"symmetricNscores1matrix\" is used, or
                  \"symmetricNscores0matrix\" when strings of Ns are detected.
 -i(gnoreindels)  Do not add  %insertion and %deletion to divergence estimate.
 -u(naligned) #   (default 35 bp or 8 AA) Maximum unaligned residues on each site of an entry for it to be removed
                  ( For proteins:(6 + (opt_d - divergence)) AA more if combined unaligned AA is smaller than that of matched protein ) 
 -z(naligned) #   (default 2x opt_unaligned) Maximum unaligned bp for both sites combined; running out of letters
 -v(naligned) #   (default off) Maximum unaligned bp on each site of the entry causing the other to be removed
                  Intended to be used when individual units and mosaic elements (e.g. FLAM and AluJ) are in input
 -ma(trix)        A path, or the name of a matrix in Matrices/crossmatch/ next to this
                  script. Default for proteins is BLOSUM62x
                  Default for DNA is \"Nscores1matrix\"  which equals \"matrix\" but does not treat N as a mismatch
                  When strings of Ns are expected (genomic assemblies, padded insertions) choose -matrix matrix
                  This is the symmetrized 25p49 matrix with N scoring -1 to everything. 
 -f(ragments)     Assumes sequences are fragments of uniquely labeled sequences. Overlapping fragments are merged by coordinate.
                  NOTE: Fragment coordinates need to be all 1-based or all half-open for merges to be perfect
 -o(verlap) #     Minimum overlap before merging fragments. Default 1 (abutting seqs stay separate). 
                  If negative, neighboring seqs are joined with N strings.
 -n(o_crossmatch) Skip the cross_match comparison; only useful with -f or for removing identical (same orientation) seqs 
 -me(rge)         Merge fragments overlapping by less than cutoff div (only necessary when -fragments is not used) 
 -s(core) #       Score to report. Defaults are 100 for proteins and 1000 for DNA (~111 bp perfect match, 130 bp 5% divergence, 160 bp 10% div)
 -w(ordlength) #  Minmatch (default 25 for DNA, 5 for AA)
 -r(eanalyze)     Takes existing selfcomparison file (<file>.vse) to rerun with different cutoff parameters.
                  Expects a modified sequence file (<file>.DADmod) if in the original run fragments were joined,
                  perfect duplicates were removed, etc. as indicated in <file>.DADmod.mods
 -j(usthalf)      Removes a sequence with matches from <= -u of either end covering half or more of the sequence    
By default the second sequence in a cross_match comparison of equals is removed. This can be modified by either
 -a(mbigfilter)   Removes sequence with higher % of ambiguous bases or AA (by default on with opt_protein)
 -t(able) <file>  Table with sequence_name,score for all sequences. The lower scoring match will be removed 
 -k(eep)          Fasta file or list (one sequence name per line) of sequences (e.g. well curated sequences) 
                   to be kept no matter, except when occurring twice with the same sequence AND name. Does not work with -merge.

Note: some of these options only work on DNA sequences.\n";

my $commandline = join " ", $0, @ARGV;

my @opts = qw (ambigfilter blast comparative=s divergence=s fragments ignoreindels justhalf keep=s matrix=s merge no_crossmatch overlap=i protein reanalyze score=s table=s znaligned=s unaligned=s vnaligned=s wordlength=i);
our ($opt_ambigfilter, $opt_blast, $opt_comparative, $opt_divergence, $opt_fragments, $opt_ignoreindels, $opt_justhalf, $opt_keep, $opt_matrix, $opt_merge, $opt_no_crossmatch, $opt_overlap, $opt_protein, $opt_reanalyze, $opt_score, $opt_table, $opt_znaligned, $opt_unaligned, $opt_vnaligned, $opt_wordlength);
unless (&GetOptions(@opts) && $ARGV[0]) {
  print STDERR "\n$USAGE";
  exit (1);
}

my $div = 10;
$div = $opt_divergence if $opt_divergence;
if ($opt_comparative && ($opt_comparative !~ /^\d\d/ || $opt_comparative < 60 || $opt_comparative > 100) ) {
  die "Choose a value between 60 and 100 for option -comparative. Good cutoff values are in the 85-95 range.\n";
}
my $unaligned = 35;
$unaligned = 8 if $opt_protein;
if (defined $opt_unaligned) {
  die "Option unaligned requires a whole positive number\n" if $opt_unaligned !~ /^\d+$/ || $opt_unaligned == 0;
  $unaligned = $opt_unaligned;
}
my $totalunaligned = 2* $unaligned;
if (defined $opt_znaligned) {
  die "Option totalunaligned requires a whole positive number\n" if $opt_znaligned !~ /^\d+$/ || $opt_znaligned == 0;
  $totalunaligned = $opt_znaligned;
}
my $vnaligned = 0;
if ($opt_vnaligned) {
  die "Option unaligned requires a whole positive number\n" if $opt_vnaligned !~ /^\d+$/;
  $vnaligned = $opt_vnaligned + 1;
}

my $merge = "";
$merge = 1 if $opt_merge;
my $overlap = 0;
$overlap = ($opt_overlap -1) if defined $opt_overlap;
# used in 1-based system;

my $file = $ARGV[0];
die "Cannot open $file\n\n$USAGE" unless -s $file;

# Check what this run needs before any output files are written.
unless ($opt_no_crossmatch) {
  if ($opt_protein) {
    &needonpath("cross_match");
  } else {
    die "Cannot find Xmatch.pl at $xmatch\n" unless -x $xmatch;
    if ($opt_blast) {
      &needonpath("rmblast.pl","rmblastn","makeblastdb");
    } else {
      &needonpath("cross_match");
    }
  }
}

my %qualscore = ();
my $tablecheck = "";
my $poundthere = ""; # the table fil may contain the name of the sequences without the class indication

my $matrix = "$matrixdir/symmetricNscores1matrix";
$matrix = "$matrixdir/BLOSUM62x" if $opt_protein;
if ($opt_matrix) {
  $matrix = -s $opt_matrix ? $opt_matrix : "$matrixdir/$opt_matrix";
  die "Cannot open $opt_matrix, nor $matrixdir/$opt_matrix\n" unless -s $matrix;
}

if ($opt_table) {
  die "Cannot read the table file $opt_table\n" unless -s $opt_table;
  open (INTABLE, $opt_table);
  while (<INTABLE>) {
    if ( /^\s*(\S+)\s*,\s*([\d\.]+)\s*/ ) {
      $qualscore{$1} = $2;  
      my $sc = $2;
      my $herename = $1;
# sometimes I use the name including #classification, sometimes just
# the name itself. A mistake, but harmles unless there is an entry
# that only differs after the #
      ++$poundthere if $herename =~ s/\#\S+//;
      ++$tablecheck;
      $qualscore{$herename} = $sc;
    }
  }
  die "$USAGE\nCould not read the format in table file $opt_table\n" unless $tablecheck;
  if ($poundthere && $tablecheck > $poundthere) {
    warn "Warning: not all entries in the table appear to have a class indicated\n";
  }
  # Check if all entries in the fasta file have occured in the opt_table file
  my $missingnames;
  my $allthefastalines = `grep \"^>\" $file`;
  while ($allthefastalines =~ s/^>(\S+)//m) { # /m makes it anchor per line
    my $fullname = $1;
    unless ( defined $qualscore{$fullname} ) {
      if ($poundthere && $tablecheck > $poundthere) {
	$fullname =~ /^(\S+)\#/;
	unless ( defined $qualscore{$1} ) {
	  $missingnames .= "$1,0\n";
	}
      } else {
	$missingnames .= "$fullname,0\n";
      }
    }
  }
  if ($missingnames) {
    warn "Some entries in $file do not have a quality score in $opt_table.\nThese are (with quality zero given for now):\n";
    die "$missingnames";
  }
}

my %keep;
if ($opt_keep) {
  die "$USAGE\nCannot read the file $opt_keep with names of entries to keep.\n" unless -e $opt_keep;
  if (-s $opt_keep) {
    open (INKEEP, $opt_keep);
    while (<INKEEP>) {
      # one can make the table by just grepping '>' from a fasta file
      if (/^>(\S+)/) { # (grepped) headers of a fasta file
	++$keep{$1};
      } elsif (/^(\S{1,50})\s/) { # just a list of names
	# does not really matter if some of these are unexpectedly short sequences; they won't match sequence names (one would hope)
	++$keep{$1};      
      } 
    }
    my $nrkeep = keys %keep;
    if ($nrkeep) {
      warn "$nrkeep sequences in $opt_keep are excluded from exclusion\n";
    } else {
      die "$USAGE\nCould not read any sequence names in the file $opt_keep\n";
    }
  }
}

my $modseqfile = "$file.DADmod"; # will contain modified sequences (empty and duplicate sequences removed; fragments joined)
my $modnotesfile = "$modseqfile.mods";
if ($opt_reanalyze) {
  my $compfile = "$file.vse";
  die "Missing $compfile required for the option -reanalyze\n" unless -s $compfile;
  if ($opt_comparative) { 
 # check if the .vse file indeed contains output of showing relative
 # per-base-score Checking the file $file.DADinfo is not reliable. The first
 # line is the previous command line, but this could already be a reanalysis in
 # which the user did not type -c (while the vse file does contain the necessary
 # info to do a -c run)
    open (IN, $compfile);
    while (<IN>) {
      if (/\(\d+\)/) {
	my @bit = split;
# expect (output of Xmatch.pl -p)
# 100.  2704  0.00  0.00  0.00  DF0282119.1   1 304 (0)    ManPen-1.181_LTR#LTR/ERVL     1  304 (0)
# or
# 99.7  3125  0.00  0.00  0.00  PerLon-1.105#LTR/ERVL   1 349 (0)    DF0279966.1               1 349 (0)
# not the default output:
# 12928  0.14  0.00  0.00  DF0275588.1#LINE/L1     1 1443    (0)    ChiLan-1.140#LINE/L1     1 1443    (0)
	unless ($bit[4] =~ /^\d+\.\d\d$/ && $bit[8] =~ /^\(\d+\)/) {
	  die "The option -comparative can only be used in reanalysis when the original run was done with that option.
The file $compfile does not seem to contain the expected output for a -c run.\n";
	}
	last;
      }
    }
  }
  if (-e $modnotesfile ) {
    if (-s $modnotesfile ) { # something was written to it
      die "Missing $modseqfile with modifications indicated in $modnotesfile (required for the option -reanalyze).\n" unless -s $modseqfile;
      open (IN, $modseqfile);
    } else {
      open (IN, $file);
    }
  } else { 
    die "A file $modnotesfile is expected if a previous DeleteAlmostDupsfromSeqfile.pl on this file was performed.\n
Without that the option reanalyze cannot work.\n"
  }
  my $prevfile = $file."\.selection";
  my $newfilename = $prevfile."_previous";
  rename ($prevfile, $newfilename) if -s $prevfile;
} else {
  open (IN, $file); # seq file
}
my ($name,$count,$lastentry,$Nstringdetected);
my @names = ();
my %entry; # values: fasta line + sequence; keys: sequence name
my (%id, %nr,%seq,%data);
my (%orderN, %orderE); # order of the sequences in the original file with Names or Entries (fastaline & sequence) as keys
my (%core,%firstnr,%corrections,%ambigcount,%length);
my %exactas; #keys: remaining entry names; values: lists of names of entries with identical sequences as key
my %mention; #keys: remaining entry names; values: lists of names of non-idential entries deleted because of the key
my (%relpairsc,%pairdiv); # hashes keeping the relative score (in $opt_comparative mode) and divergence level of pairs of seqs
while (<IN>) {
  if (/^>(\S+)/) {
    $name = $1;
    # At the line "if ($opt_fragments && $lastname =~ /^$corename[_:]\d/) {"
    # names (at tha point $corename) containing a backslash die with the following warning
    # " Character following \P must be '{' or a single-character Unicode property name in regex; marked by <-- HERE in m/^Damb\P- <-- HERE element_T[_:]\d/ at /home/asmit/bin/DeleteAlmostDupsfromSeqfile.pl line 503, <IN> line 15294.
    # quotemeta doesn't help; there must be a better solution
    die "$name contains a backslash, which the script currently cannot handle. Please remove/replace all backslashes.\n" if $name =~ /\\/;
    $orderE{$lastentry} = $count if $count;
    ++$count;
    my $naam = $name;
    if ($nr{$name}) {
      $name .= "..$nr{$naam}";
      $name =~ s/(\#\S+)(\.\.\d+)$/$2$1/;
      # an unlikely string for a name or accession number
      while ($entry{$name}) {
	++$nr{$naam};
	$name .= "..$nr{$naam}";
	$name =~ s/(\#\S+)(\.\.\d+)$/$2$1/;
      }
      s/>$naam/>$name \($naam\)/;
    } else {
      push @names, $name;
    }
    # If reading in a file with previously deredundanted entries
    # Common with $opt_reanalyze reading in the modseqfile
    if ( s/\(a.k.a.\s+(\S+)\s*\)// ) {
      $exactas{$name} = "$1,";
    }
    if ( s/Kept in favor of\s+(\S+)// ) {
      $mention{$name} = "$1,";
    }
    if ($opt_table && $name =~ /\.\.\d+(\#\S+)?$/) {
      my $oldname = $name;
      unless ($oldname =~ s/\.\.\d+\#/\#/) {
	$oldname =~ s/\.\.\d+$//;
      }
      $qualscore{$name} = $qualscore{$oldname};
    }
    /^>\S+(.*)\n/; # . does not include \n
    my $stuff = $1;
    if ($stuff && $stuff =~ /\S/) {
      $data{$name} = $stuff;
      $data{$name} =~ s/^\s+//;
      $data{$name} =~ s/\s+$//;
    }
    ++$nr{$naam}; 
    $ambigcount{$name} = 0; # to be set in sequence
    $id{$name} = $name;  #looks weird, but changes below
    $orderN{$name} = $count;
    $corrections{$name} = 0;
    if ($opt_protein && ($_ =~ /frameshift/ || $_ =~ /stopcodon/) ) {
      if ( /(\d+)\s+frameshift/ ) {
	$corrections{$name} += $1;
      }
      if ( /(\d+)\s+stopcodon/ ) {
	$corrections{$name} += $1;
      }
    }
  } else {
    chomp;
    # change to upper case so that entries are not differentiated solely by case
    tr/a-z/A-Z/;
    $seq{$name} .= $_;

    # count ambiguous AA or bp in sequence
    if ($opt_protein) {
      $ambigcount{$name} += tr/xX//;
    } else {
      $ambigcount{$name} += (tr/nrykmswNRYKMSW//); 
    }
    $length{$name} += (length $_);
  }
  $entry{$name} .= $_;

  #override the default matrix when long strings of Ns occur in the query
  unless ($matrix !~ /Nscores1matrix$/ || $Nstringdetected) {
    ++$Nstringdetected if $entry{$name} =~ /[Nn]{50,}/;
  }
  $lastentry = $entry{$name};
  if ( $name =~ /^(\S+)[:_](\d+)[-_]\d+(_[\dRC])?(_\d)?$/ ) {
    $core{$name} = $1;
    $firstnr{$name} = $2;
  } else {
    $core{$name} = $entry{$name};
    $firstnr{$name} = 0;
  }
}
$orderE{$lastentry} = $count if $count;

if ($Nstringdetected) {
  $matrix = "$matrixdir/symmetricNscores0matrix";
  warn "Strings of Ns longer than 50 detected in query. 
Changing the matrix to $matrix to avoid false matches.\n"
}

my %ambigfrac;
if ($opt_ambigfilter) {
  foreach my $sqname (keys %length) {
    if ($length{$sqname}) {
      $ambigcount{$sqname} = 0 unless $ambigcount{$sqname};
      $ambigfrac{$sqname} =  $ambigcount{$sqname}/$length{$sqname};
    } else {
      warn "Sequence $sqname has no length";
    }
  }
}


unless ($opt_reanalyze) {
  open (OUTNOTES, ">$modnotesfile");
  # this file created so that when option -r is used in a possible
  # next attempt with different parameters the out from this first
  # cleanup is not lost

  # remove exact duplicates (i.e. name and sequence are the same)
  # and rename those with same name but different sequence
  foreach my $naam (@names) {
    if ($nr{$naam} > 1) { 
      for (my $i = 1; $i < $nr{$naam}; $i++) {
	my $nom = "$naam..$i";
	$nom =~ s/(\#\S+)(\.\.$i)/$2$1/;
	my $j = $i -1;
	if ($seq{$nom} =~ /[ACGTacgt]/) {
	  if ($seq{$nom} eq $seq{$naam}) {
	    print OUTNOTES "Removed a duplicate of entry $naam\n";
	    delete $entry{$nom};
	  } else {
	    while ($j > 0) {
	      my $nome = "$naam..$j";
	      $nome =~ s/(\#\S+)(\.\.$j)/$2$1/;
	      if ($entry{$nome} && $seq{$nom} eq $seq{$nome}) {
		print OUTNOTES "Removed a duplicate of entry $nome\n";
		delete $entry{$nom};
		last;
	      }
	      --$j;
	    }
	    if ($entry{$nom}) {
	      print OUTNOTES "Renamed an entry with the previously used name $naam for a different sequence $nom\n";
	      
	    }
	  }
	} else {
	  print OUTNOTES "Removed empty sequence $nom\n";
	  delete $entry{$nom};
	}
      }
    }
  }

  my %dup;
  my $lastname = "";
  my $lastbeg = 0;
  my $lastend = 0;
  my $lastrev = "";
  my @entries;
  my @sortedentries = sort { $core{$a} cmp $core{$b} || $firstnr{$a} <=> $firstnr{$b} } keys %entry;
  foreach my $naam (@sortedentries) {
    my ($corename, $begin, $eind) = (0,0,0);
    my $rev = "F";
    my $seq = $entry{$naam};
    $seq =~ s/^>(\S+).*\n//;
    if ($dup{$seq}) {
      # the very same sequence has already been observed; 
      # $dup{$seq} contains the name ($naam) of that previous entry
      my ($nm1,$nm2) = ($naam,$dup{$seq});
      $nm1 =~ s/\#\S+// unless defined $qualscore{$naam};
      $nm2 =~ s/\#\S+// unless defined $qualscore{$dup{$seq}};
      if ( ( !$keep{$dup{$seq}} || $keep{$naam} ) && 
	  ( $opt_table && $qualscore{$nm1} > $qualscore{$nm2} || 
	    !$keep{$dup{$seq}} && $keep{$naam} ) ) {
	my $despite = "";
	$despite = "despite both names being on the keeplist $opt_keep, as $naam had a higher score in $opt_table" if $keep{$dup{$seq}};
	my $line = "Removed $dup{$seq} as an exact duplicate of $naam $despite\n";
	print OUTNOTES $line;
	delete $entry{$dup{$seq}};
	$exactas{$naam} .= "$dup{$seq},";
	$exactas{$naam} .= "$exactas{$dup{$seq}}," if $exactas{$dup{$seq}};
	# later on: $dup{$seq} = $naam;
      } else {
	# the second entry encountered is removed if there is no keeplist nor a
	# score table, no keeplist and they have the same score in table, or
	# both on keeplist and (no table or 2nd table score <= 1st)
        my $despite = "";
        $despite = "despite both names being on the keeplist $opt_keep" if $keep{$name};
	my $line = "Removed $naam as an exact duplicate of $dup{$seq} $despite\n";
        print OUTNOTES $line;
	delete $entry{$naam};
	$exactas{$dup{$seq}}.= "$naam,";
	$exactas{$dup{$seq}}.= "$exactas{$naam}," if $exactas{$naam};
	next;
      } 
    }
    if ($opt_fragments) {
      $corename = $naam;
      # SplitSeqs used to add _begin_end or _begin_end_R to the name duplicate
      # entries get a _nr prefix from above. These shouldn't get here,
      # but just to be sure a _(single digit) suffix is included
      # fortunately single digit ends don't happen in fragments
      #Example:
      #>fullTreeAnc238refChr12149_10771_11076_C
      if ($corename =~ s/_(\d+)_(\d+)(_[\dRC])?(_\d)?$//) { 
	$begin = $1;
	$eind = $2;
	if ($3 && $3 =~ /[RC]/) {
	  $rev = 'R';
	} elsif ($begin > $eind) {
	  $rev = 'R';
	  my $temp = $begin;
	  $eind = $begin;
	  $begin = $temp;
	}
	# We now use the format name:begin-end_R 
      } elsif ($corename =~ s/:(\d+)\-(\d+)(_[\dRC])?(_\d)?$//) {
	$begin = $1;
        $eind = $2;
	if ($3 && $3 =~ /[RC]/){
          $rev = 'R';
	} elsif($begin> $eind) {
          $rev = 'R';
          my $temp = $begin;
          $eind = $begin;
          $begin = $temp;
        }
      } else {
	# if a full sequence has been collected SplitSeqs now removes begin and end info
	# we're hacking it by providing a begin of 1 and an extreme end
	$begin = 1;
	$eind = 999999999;
	$rev = 'R' if $corename =~ /_[RC](_\d)?$/;
	print OUTNOTES "Added fictional begin and and to the name of the full-length contig $corename\n";
      }
      unless ($begin =~ /^\d+$/) {
	die "Expected a number (beginning position) for $begin $naam $corename $begin $eind $rev\n";
      }
    }
    if ($opt_fragments && $lastname =~ /^$corename[_:]\d/) {
      # two fragments of the same database entry
      if ($lastbeg <= $begin && $lastend >= $eind ) {
	print OUTNOTES "Removed $naam because of $lastname\n";
	$entry{$naam} = "";
	$begin = $lastbeg;
	$eind = $lastend;
	$naam = $lastname;
      } elsif ( $lastbeg >= $begin && $lastend <= $eind ) {
	pop @entries;
	$entry{$lastname} = "";
	print OUTNOTES "Removed $lastname because of $naam\n";
	push @entries, $entry{$naam};
	$dup{$seq} = $naam;
      } elsif ( $begin <= ($lastend-$overlap) && $eind >= ($lastbeg + $overlap) ) {
	# note: by default $overlap = 0 and flush fragments will not be merged 
	my $newname = $naam;
	my $strictoverlap = $overlap + 100; 
	# strict overlap as we often do not want to merge barely/non-overlapping frags in opposite orientation 
	if ($rev ne $lastrev && $begin <= ($lastend-$strictoverlap) && $eind >= ($lastbeg + $strictoverlap) ) {
	  if ($lastrev eq 'F') {
            $entry{$naam}  =~ /^(>.+)\n(\w+)/;
            my $newfaline = $1;
            my $revseq = reverse $2;
            $newfaline =~ s/^(\S+)_[CR](_\d)?/$1/;
	    # $newname will be used to build the new fasta line
	    $newname =~ s/^(\S+)_[CR](_\d)?/$1/;
            $revseq =~ tr/ACGTYRMKSWHBVD/TGCARYKMSWDVBH/;
            $revseq =~ tr/acgtyrmkswhbvd/tgcarykmswdvbh/;
            $entry{$naam}  = "$newfaline\n$revseq\n";
            $rev = 'F';
	  } elsif ($rev eq 'F') {
	    $entry{$lastname}  =~ /^(>.+)\n(\w+)/;
	    my $newfaline = $1;
	    my $revseq = reverse $2;
	    $newfaline =~ s/^(\S+)_[CR](_\d)?/$1/;
	    $revseq =~ tr/ACGTYRMKSWHBVD/TGCARYKMSWDVBH/;
	    $revseq =~ tr/acgtyrmkswhbvd/tgcarykmswdvbh/; 	    
	    $entry{$lastname}  = "$newfaline\n$revseq\n";
	    $lastrev = 'F';
	  } else {
	    die "\$rev ($rev) ne \$lastrev ($lastrev) but neither eq F\n";
	  }
	}
	if ($rev eq $lastrev) {
	  if ($lastbeg < $begin) {
            # entries were once sorted alphabetically so that "later"
            # fragments sometimes were sorted before earlier,
            # e.g. seq_1000_2000 comes before seq_900_2000. This is
            # not the case anymore so $lastend should never be > $eind
            # For code to handle that situation see
            # ~/bin/Olderversions/DeleteAlmostDupsfromSeqfile_210831.pl
	    pop @entries;
	    my $data = &getdata($data{$lastname},$data{$naam});
	    my $Ns = "";
	    my $dist = $begin - 1 - $lastend;
	    # begin 1-based
	    $Ns = "N"x $dist if $dist > 0;
	    if ($rev eq 'F') {
	      # --->   or --->
	      #  --->         --->
	      print OUTNOTES "1 $lastbeg $begin Merged $lastname and $naam\n";
	      $newname =~ s/([_:])$begin([_-])/$1$lastbeg$2/;
	      $entry{$lastname}  =~ s/^>.+\n//;
              $entry{$naam}  =~ s/^>.+\n//;
	      if ($dist >= 0) {
		$entry{$newname} = ">$newname$data\n$entry{$lastname}$Ns$entry{$naam}";
	      } else {
		my $begdif = $begin - $lastbeg;
		my $add = substr($entry{$lastname},0,$begdif);
		$entry{$newname} = ">$newname$data\n$add$entry{$naam}";
	      }
	    } else {
	      # <---  or <---
	      #  <---        <---
              print OUTNOTES "4 $begin $lastbeg Merged $lastname and $naam\n";
	      die "Oops. \$naam ($naam) is not identical to \$newname ($newname)\n" if $newname ne $naam;
	      $entry{$lastname}  =~ s/^>.+\n//;
	      $entry{$naam} =~ s/^(>\S+[_:])$begin([_-])/$1$lastbeg$2/;
	      $entry{$naam} =~ s/^(>\S+).*/$1$data/; #.* and not .+
	      chomp $entry{$naam};
	      # I think $newname is always identical to $naam here
	      # (only different when the orientation of the last entry
	      # was changed from _R to forward to make the merge possible
	      if ($dist >= 0) { # non-overlapping
		$entry{$naam} .= "$Ns$entry{$lastname}";
	      } else {
## $begdif used to have 1 added to it, for the end of line in the flipped sequence.
## That no longer seems necessary; perhaps the lastname entry was chomped, or
## there was a mix of half-open and 1-based entries when I first tested it.
		my $begdif =  $begin - $lastbeg;
		my $add = substr($entry{$lastname},-$begdif); 
		$entry{$naam} .= "$add"; # $add includes the line break
	      }
## the following seems necessary, even if I didn't have that n the old code
	      $entry{$lastname} = "";
	      
	    }
	    $begin = $lastbeg;
	  } else {
	    warn "How did $lastname start after $naam?\n";
	  }
	  $orderE{$entry{$newname}}= $orderN{$lastname}; 
# the entry of lastname has been messed up, but the numbering is the same for the names as the entries
	  $orderN{$newname} = $orderN{$lastname}; # in case the next entry also gets merged;
	}
	if ($newname ne $naam) {
	  push @entries, $entry{$newname};
	  $entry{$lastname} = "";
	  $entry{$naam} = "";
	  $naam = $newname; # to create proper $lastname
	} else {
	  push @entries, $entry{$naam};
	}
	$dup{$seq} = $naam;
      } else {
	push @entries, $entry{$naam};
	$dup{$seq} = $naam;
      }
    } else {  
      push @entries, $entry{$naam};
      $dup{$seq} = $naam;
    }
    if ($opt_fragments) {
      $lastname = $naam;
      $lastbeg = $begin;
      $lastend = $eind;
      $lastrev = $rev;
    }
  }
  %dup = ();
  close OUTNOTES;
  if (-s $modnotesfile) {
    open (OUT, ">$modseqfile");
    foreach my $nom (sort { $orderN{$a} <=> $orderN{$b} } keys %entry) {
      if ( $exactas{$nom} ) {
	$exactas{$nom} =~ s/,$//;
	$exactas{$nom} =~ s/,,/,/g;
	$entry{$nom} =~ s/^(>.+?)\n/$1 \(a.k.a. $exactas{$nom}\)\n/;
      }
      $entry{$nom} .= "\n" unless $entry{$nom} =~ /\n$/;
      print OUT "$entry{$nom}";
    }
    close OUT;
  }
}
my $infofile = "$file.DADinfo";
if (-s $infofile) {
  warn "Renaming existing file $infofile $infofile.before\n";
  rename $infofile, "$infofile.before";
}
open (OUTINF, ">$infofile");
print OUTINF "$commandline\n\n";
if (-s $modnotesfile) {
  my $notes = `cat $modnotesfile`;
  print OUTINF "$notes\n";
  if ($opt_no_crossmatch) {
    print OUTINF "With option -n(o_crossmatch), no sequence alignments were performed.\n";
  } elsif ($opt_blast) {
    print OUTINF "Following self comparison with rmblastn:\n\n";
  } else {
    print OUTINF "Following self comparison with cross_match:\n\n";
  }
}

my $minscore = 1000;
if ($opt_score) {
  $minscore = $opt_score;
} elsif ($opt_protein) {
  $minscore = 100;
}
my $query = $file;
$query = $modseqfile if -s $modnotesfile;
unless ($opt_no_crossmatch) {
  my $tempvse = "$file.vse";
  if ($opt_reanalyze) { # the .vse file already exists
    warn "Renalyzing $tempvse\n";
  } else { # make the .vse file
    my $sortedvse = "$tempvse".".sorted";
    open (OUTSORTED, ">$sortedvse");
    # For partially masked sequences it is key to use one of my own
    # matrices, as Ns or Xs are otherwise scored negative and a high
    # mismatch level is returned
    if ($opt_protein) {
      my $minmatch = 5;
      $minmatch = $opt_wordlength if $opt_wordlength;
      my $comd = "cross_match $query -masklevel 101 -minmatch $minmatch -minscore 100 -matrix $matrix -gap_init -12 -gap_ext -2 -alignments";
      warn "\n$comd\n";
      print OUTSORTED "$comd\n\n";
      system "$comd > $tempvse";
    } else {
      my $minmatch = 25;
      $minmatch= $opt_wordlength if $opt_wordlength;
      my $comd = "$xmatch $query -level 101 -minmatch $minmatch -score $minscore -matrix $matrix -gap_init -30 -extension -6 -bandwidth 100 -q ";
      if ($opt_comparative) {
	# overrides the $opt_matrix choice
	$comd .= " -perbasescore -cg";
      } elsif ($opt_blast) {
	$comd .= " -blast";
      }
      warn "\n$comd\n";
      print OUTSORTED "$comd\n\n";
      system "$comd > $tempvse";
    }
  
# clean up cross_match output and sort by start position, score and divergence
# so that the best match tends to end on top (except when it starts at a later base ...)
    if ($opt_comparative) {
      warn "Sorting cross_match-formatted output by score relative to self-score\n";
    } else {
      warn "Sorting cross_match-formatted output by score per length of query\n";
    }
    my (@xmlines,%linescore,%linediv,%nm,%linebeg,%aligntext,%relscore);
    my ($skipthisjunk,$cmline,$lastquery);
    open (XM, $tempvse); #crossmatchfile comp against self, masklevel 101 
    while (<XM>) {
      if ( /^Num. pairs/ || $_ =~ /^Discrepancy summary/) {
	$skipthisjunk = 1;
      } elsif ( $skipthisjunk && $_ =~ /\(\d+\)/ ) {
	$skipthisjunk = 0;
	print OUTSORTED "\n";
      }
      next if $skipthisjunk;
      if (/^Exact duplicate|\(perfect\)/) {
	# cross_match also reports perfect reverse complemented matches and does
	# not use the second sequence in the "assembly"
	print OUTSORTED;
      } elsif ( /\(\d+\)/ ) {
	my @bit = split;
	next unless $bit[9];
	if ( $bit[8] =~ /^\(\d+\)$/ ) { # with $opt_comparative the first number is an extra relative score
	  $relscore{$_} = shift @bit;
	} else {
	  if ($opt_comparative) {
	    if ($opt_reanalyze) {
	      warn "The file $tempvse does not contain all information needed for the -comparative mode in this reanalysis\n";
	      $opt_comparative = "";
	    } else {
	      die "Unexpected cross_match line for -comparative mode:\n$_";
	    }
	  }
	  $bit[7] =~ tr/()//d;
	  # create a poor man's version of the relative score to sort alignments by
	  $relscore{$_} = $bit[0]/($bit[6]+$bit[7]);
	}
	next if $bit[0] < $minscore; # can happen when reanalysis is done with a higher minscore
	if ($bit[8] eq 'C' && $bit[4] eq $bit[9] || $bit[11] =~ /^\(\d+\)$/ && $bit[4] eq $bit[8] ) {
	  next; # a match against self. Otherwise many hairpin elements would be deleted as well as short tandem repeat models.
	}
# changed all sorting to the way opt_comparative lines are sorted
# which has the benefit of seeing the closest matches first
	$linescore{$_} = $bit[0];
	$linediv{$_} = $bit[1];
	$nm{$_} = $bit[4];
	$linebeg{$_} = $bit[5];
	push @xmlines, $_;
	$cmline = $_;
	$lastquery = $bit[4];
      }
      $aligntext{$cmline} .= $_ if $cmline;
    }
    my (%linecount,%lastscore);
    my @sortedlines = sort { $relscore{$b} <=> $relscore{$a} ||
				 $linescore{$a} <=> $linescore{$b} || 
# the thought is, given a tie, removing small fragments first might be better
				 $linebeg{$a} <=> $linebeg{$b} ||
				 $linediv{$a} <=> $linediv{$b} } @xmlines;
    @xmlines = ();
    foreach my $xmline (@sortedlines) {
      ++$linecount{$nm{$xmline}};
      if ($linecount{$nm{$xmline}} <= 12 || # keeps the file size down; how likely is it that a query would be deletd based on the 13+ best match?
	  $linescore{$xmline} >= 0.99*$lastscore{$nm{$xmline}} ) { # 
	print OUTSORTED "$aligntext{$xmline}";
	$lastscore{$nm{$xmline}} = $linescore{$xmline} if $linecount{$nm{$xmline}} <= 12;
      } else {
	# the cross_match lines are needed though for checks like sub checkifokay (did a previously removed seq match the current bully)
	$aligntext{$xmline} =~ /^(\s*\d.+?)\n/;
	print OUTSORTED "$1\n\n";
      }
    }
    (@xmlines,@sortedlines,%aligntext) = (); # I supposed these get flushed out of memory anyway when leaving this loop
    close OUTSORTED;
    close XM;
    rename $sortedvse,$tempvse
  } # unless $opt_reanalyze

#keys: "query name"_"targetname" of aligned pair. Values: alignment score/maximum possible score for query
  my %totscore; #keys: sequence names, values: total score of all alignments with this sequence as query
  my $checksubfrags;
  open (XM, $tempvse);
  while (<XM>) {
    if (/^\s*\d/) {
      my @bit = split;
      next unless $bit[10];
      my $relsc; # relative score or SW score
      $relsc = shift @bit if $bit[4] =~ /^\d+\.\d+$/ && $bit[8] =~ /\(\d+\)/;
      # in principal one should be able to do a reanalysis not using the
      # -comparative option even if the first run was done so.
      # However, elsewhere in the code, I assume that a reanalysis is done in
      # the same mode as the original. I need to fix that.
      next if $bit[0] < $minscore; # in case reanalysis was done with a higher cutoff score
      my $pair = "$bit[4]"."_$bit[8]";
      $pair = "$bit[4]"."_$bit[9]" if $bit[8] eq 'C' || $bit[8] eq '+';
      if ($opt_comparative) {
	$relpairsc{$pair} = $relsc unless $relpairsc{$pair} && $relpairsc{$pair} > $relsc;
	$totscore{$bit[4]} += $bit[0];
      } else {
	$bit[7] =~ tr/()//d;
	if ($bit[5] - 1 <= $unaligned && $bit[7] <= $unaligned
	    && $bit[5] + $bit[7] - 1 <= $totalunaligned ) {
	  $pairdiv{$pair} = $bit[1];
	} else {
	  $pairdiv{$pair} = $div + 10;
	}
	$totscore{$bit[4]} += $bit[0];
      }
      # for long alignments we want to be sure that no subregions show high
      # divergence for example, LINE can look very similar over the conserved
      # ORF2 region bu can differ a lot in the 3' UTR and 5' region
      if ($bit[6] - $bit[5] > 1500 &&
	  ($opt_comparative && $relsc >= $opt_comparative ||
	   !$opt_comparative && $pairdiv{$pair} <= $div)) {
	# one can also skip the check if the overall divergence level is so low
	# or relative score is so high that, given all mismatches are in one
	# fragment that fragment still would pass the test.
	$checksubfrags = 1;
      }	 else {
	$checksubfrags = 0;
      }
      # currently nothing is happening but the intention is to use the existing
      # alignments in this and calculate the score for fractions of it
      # not sure if without the -c option there will always be alignments
    }
  }

  open (XM, $tempvse); #crossmatchfile comp against self, masklevel 101
  LOOP: while (<XM>) {
    if (/\(perfect\)\n/) {
      my @fields = split;
      my $nr1 = $fields[0];
      my $nr2 = $fields[1];
      $nr2 = $fields[2] if $nr2 eq 'C';
      if ( $keep{$nr2} ) {
	if ( $keep{$nr1} ) {
	  print OUTINF "Identical entries $nr1 and $nr2 retained as they are both in the \"keep these\" file $opt_keep\n";
	} else {
	  # you'll only find $nr1 in the cross_match output 
	  # so treat the entry (fasta line and sequence) of $nr2 as if it is $nr1 
	  $entry{$nr1} = $entry{$nr2};
	  # When $fields[1] eq 'C' the entries are reverse complements
	  # I believe this has no effect on the analysis yet (alignments are not
	  # inspected for orientation), but keep this in mind
	  $id{$nr1} = $id{$nr2}; # a way to scratch the first sequence name from the list 
	  chomp;
	  print OUTINF "1goneK $_ ($nr2 in keeplist)\n";
	  $exactas{$nr1} .= "$nr1,"; # counterintuitive, but identity of nr1 will be switched to nr2 later 
	  $entry{$nr2} = "";
	}
      } elsif ( $keep{$nr1} ) {
	chomp;
	print OUTINF "2goneK $_ ($nr1 in keeplist)\n";
	$exactas{$nr1} .= "$nr2,";
	$entry{$nr2} = "";
      } elsif ( $opt_table) {
	# the following shouldn't happen anymore, given the 
	# "Check if all entries in the fasta file have occured in the opt_table file " above
	if ( $qualscore{$nr2} <= $qualscore{$nr1} ) {
	  chomp;
          print OUTINF "2goneT $_ ($nr1 scores higher in $opt_table)\n";
	  $exactas{$nr1} .= "$nr2,";
	} else {
	  $entry{$nr1} = $entry{$nr2};
	  $id{$nr1} = $id{$nr2};
	  chomp;
          print OUTINF "1goneT $_ ($nr2 scores higher in $opt_table)\n";
          $exactas{$nr1} .= "$nr1,";
	}
	$entry{$nr2} = "";
      # for sequences downloaded from NCBI:
      } elsif ($opt_protein && $nr1 !~ /\|ref\|/  ##reference sequence 
	       && ($nr2 =~ /\|ref\|/ || 
		   $nr1 !~ /\|\w{4}_\w{5}$/  ##SwissProt entry 
		   && $nr2 =~ /\|\w{4}_\w{5}$/ ) ) {
	$entry{$nr1} = $entry{$nr2};
	$id{$nr1} = $id{$nr2};
	print OUTINF "1goneR $_";
      } else {
	print OUTINF "2goneR $_";
	$exactas{$nr1} .= "$nr2,";
      }
      $entry{$nr2} = "";
    }
    elsif (/\(\d+\)/) {
      next unless /\d+\s+\d+\s+\(\d+\)/; # At the end of a cross_match file are now
      # sometimes one or two lines with /\(\d+\)/ that are not what we want
      my ($unalileft2, $unaliright2, $firstname, $secname,
	  $nopound1name,$nopound2name,$match1frac,$match2frac);
      my @bit = split;
      unless ($bit[9] && $bit[11]) {
	warn "What is this line then?\n$_";
	next;
      }
      my $relscore = shift @bit if $bit[8] =~ /^\(\d+\)$/; #relative score in first column
      next if $bit[0] < $minscore; # in case reanalysis was done with a higher cutoff score       
      $firstname = $bit[4];
      if ($bit[8] eq 'C') {# && $bit[10] =~ /\(/ ) {
	$bit[10] =~ tr/[()]//d;
	$unalileft2 = $bit[10];
	$unaliright2 = $bit[12]-1;
	$secname = $bit[9];
	$match2frac = ($bit[11]-$unaliright2)/($bit[11]+$bit[10]);
      } else {
	splice @bit,8,1 if $bit[8] eq '+'; # GrepCrossmatch (or less likely RepeatMasker) output
	$secname = $bit[8];
	$unalileft2 = $bit[9]-1;
	$bit[11] =~ tr/[()]//d;
	$unaliright2 = $bit[11];
	$match2frac = ($bit[10]-$unalileft2)/($bit[10]+$bit[11]);
      }
      my $begUn = $bit[5] - 1; # unaligned bases before match
      my $endUn = $bit[7];     # unaligned bases after match
      $endUn =~ tr/[()]//d;
      my $subfrag = 0;
      if ($vnaligned) {
	if ( $unalileft2 >= $vnaligned || $unaliright2 >= $vnaligned ) {
	  $subfrag = "Kept1 (-v = $vnaligned) $_";
	} elsif ( $begUn >= $vnaligned || $endUn >= $vnaligned ) {
	  # bit[5] 1-based
	  $subfrag = "Kept2 (-v = $vnaligned) $_";
	}
      }
      # for DNA fragments of scaffolds with seqs separated by Ns -> often overlap fragments are created by Fragmentto ... in bin
      # for protein: had to label files manually before merging
      $nopound1name = $firstname;
      $nopound1name =~ s/\#.+$// unless $qualscore{$firstname};
      $nopound2name = $secname;
      $nopound2name =~ s/\#.+$// unless $qualscore{$secname};
      $match1frac = ($bit[6] - $begUn)/($bit[6] + $endUn);
      #if one of entries is gone, no use deleting the other entry
      $bit[1] = 100*$bit[1]/(100-$bit[3]); # cross_match reports mismatches/(end-begin+1) rather than over the aligned bases
      # this way you get huge underestimates of divergence when query has a large insertion compared to target
      my $fragsoverlap = 0;
      if ($opt_fragments) {
	my $name1 = $firstname;
	$name1 =~ s/[:_](\d+)_(\d+)(_[\dRC])?$//;
	my $beg1 = $1;
	my $end1 = $2;
	my $name2 = $secname;
	$name2 =~ s/[:_](\d+)_(\d+)(_[\dRC])?$//;
	my $beg2 = $1;
	my $end2 = $2;
# No harm if there are no fragments; no fragment overlap will be reported
	if ($beg1 && $beg2) {
	  if ($name1 eq $name2 && $bit[8] ne 'C' &&
	      ($beg1 < $end2 && $end1 > $beg2 ||
	       $beg2 < $end1 && $end2 < $beg1) ) {
	    $fragsoverlap = 1;
	  }
	}
      }
      my $divergence = $bit[1];
      $divergence += ($bit[2] + $bit[3]) unless $opt_ignoreindels;
      my $goahead = 0;
      # both entries still have to be there:
      if ( $entry{$firstname} && $entry{$secname} && 
	   # one of the entries has to be deletable:
	   (!$keep{$firstname} || !$keep{$secname} ) ) {
	if ($opt_comparative) {
	  die "No relative score for $_\n" unless $relscore;
	  if ($relscore >= $opt_comparative) {
	    $goahead = 1;
	  } else {
	    last LOOP; # as $tempvse is ordered by $relscore, forthcoming lines will all fail
	  }
	} else {
	  $goahead = 1 if $divergence < $div;
	}
      }
      if ($goahead) {
	# check if you do not remove a sequence A for which one or ore other
	# sequences were removed (stored inthe hash %mention) that are not
	# similar enough to sequence B, for which you are about to remove A.
	my $oppositepair = "$secname"."_$firstname";
	$relpairsc{$oppositepair} = 1 if $opt_comparative && !$relpairsc{$oppositepair};

	my ($okaytoremove1,$okaytoremove2) = (1,1);
	if ($keep{$firstname}) {
	  $okaytoremove1 = 0;
	} elsif ($mention{$firstname}) {
	  $okaytoremove1 = &checkifokay($firstname,$secname);
	}
	if ($keep{$secname} || 
	    $opt_comparative && $relpairsc{$oppositepair} < ( $opt_comparative - 20 ) ) { # could make this much closer to cutoff probably
	  $okaytoremove2 = 0;
	} elsif ($mention{$secname}) {
	  $okaytoremove2 = &checkifokay($secname,$firstname);
	}
	next unless $okaytoremove1 || $okaytoremove2;
	if ($opt_protein) {
	  # always forward orientation...
	  my $len1 = $bit[6] + $endUn;
	  my $len2 = $bit[10] + $bit[11];
	  my $diverdiff = $div - $divergence; 
# going to use this in a completely arbitrary way by allowing 
# more similar peptides to be merged or deleted even when they stick out a bit further
	  if ($merge && ($firstname eq $secname || $divergence < $div/5 ) && $bit[8] ne 'C') {
	    if ($bit[9] == 1 && $endUn == 0) {
	      # not sure anymore why $seq{$secname} needs to be revised
	      # it could be because it is already itself a merge product
	      ($seq{$secname} = $entry{$secname}) =~ s/^>.*\n([A-Za-z]+)/$1/;
	      print OUTINF "2mergedwith1 $_";
	      $entry{$firstname} .= substr($seq{$secname},$bit[10]);
	      $entry{$firstname} =~ s/^(>\S+)/$1 Nterm_extended_with_$bit[11]\_AA_of_$secname/;
	      $entry{$secname} = "";
	    } elsif ($begUn == 0 && $bit[11] == 0) {
	      ($seq{$firstname} = $entry{$firstname}) =~ s/^>.*\n([A-Za-z]+)/$1/;
	      print OUTINF "1mergedwith2 $_";
	      $entry{$secname} .= substr($seq{$firstname},$bit[6]);
	      $entry{$secname} =~ s/^(>\S+)/$1 Nterm_extended_with_$endUn\_AA_of_$firstname/;
	      $entry{$firstname} = "";
	    }
	  } elsif ($entry{$firstname}  && $entry{$secname}) {
# another hack: we much prefer clean translations over ones with
# ambiguities, stopcodons, or frameshifts (all presented by Xs) since
# we want to choose the longest sequence as well (often entries have
# internal deletions or ill-conceived "introns", we could subtract
# ambiguous etc AA from the length of the protein, to make the less
# ambguous etc preferable. The problem is in guessing a right weight
	    my $cor1 = $corrections{$firstname};
	    my $cor2 = $corrections{$secname};
	    my $amb1 = $ambigcount{$firstname};
	    my $amb2 = $ambigcount{$secname};
	    $len1 -= (4*$cor1 + 2.5*$amb1);
	    $len2 -= (4*$cor2 + 2.5*$amb2);
	    --$bit[9]; #make zero-based;
	    my $relaxedUn = $unaligned + 6 + $diverdiff;
# delete the "longer" of two sequences if the "shorter" one is forced to be kept by $keep{}
	    if (!$keep{$firstname} && $len2 >= $len1 && $okaytoremove1 &&
		( $begUn <= $unaligned && $endUn <= $unaligned ||
		  $begUn <= $relaxedUn && $endUn <= $relaxedUn 
		  && ($begUn + $endUn) < ($bit[9] + $bit[11]) ) ) {
	      &removeone($secname,$firstname,"1gone $cor1 $cor2 $amb1 $amb2 $_",$subfrag);
	    } elsif ( !$keep{$secname} && $len1 >= $len2 && $okaytoremove2 &&
		      ($bit[9] <= $unaligned && $bit[11] <= $unaligned ||
		       $bit[9] <= $relaxedUn && $bit[11] <=$relaxedUn 
		       && ($begUn + $endUn) >= ($bit[9] + $bit[11]) ) ) {
	      &removeone($firstname,$secname,"2gone $cor1 $cor2 $amb1 $amb2 $_",$subfrag);
	    }
	  }
	} else {  # batch of DNA sequences
	  my $remove2isfine;
	  my $score2over1 = 0;
	  if ($okaytoremove2) { # checked that any sequences removed for it align good enough to query as well
	    $remove2isfine = 1;
	    # do not remove sequence 2 (target) if its alignment against query
	    # scores < $opt_comparative % of a self match and score < 0.98 x the
	    # comparative score of the query against target
	    my $compratio = 1;
	    $compratio = $relpairsc{$oppositepair} / $relscore if $opt_comparative;
	    # let the total SW score decide, given the ratio in comparative scores
	    $score2over1 = 1 if $totscore{$secname}/$totscore{$firstname} > $compratio**2;
	    if ( $opt_comparative ) {
	      if ($relpairsc{$oppositepair} < $opt_comparative && $compratio < 0.98 
		  # if 1 vs 2 scores just above, and 2 vs 1 just below opt_comparative, there is no reason to discriminate
		  || $keep{$secname} ) {
		# can add more reasons here as well
		$remove2isfine = 0;
	      }
	    }
	  }
	  next unless $okaytoremove1 || $remove2isfine;

	  # Are the unaligned ends of firstname at least as long as those of secname - 20 bp?
	  # Are the unaligned ends of secname shorter than $unaligned?
	  # ($unaligned is set by $opt_unaligned or is 35 by default for DNA)
	  # If so secname is contained in firstname seq.
	  # Skipped that 
	  if ( $unalileft2 <= $unaligned && $unaliright2 <= $unaligned
	      && $unalileft2 + $unaliright2 <= $totalunaligned ) { # &&
	       # Are the unaligned ends of firstname at least as long as those of secname - 20 bp? 
	       # If so secname is contained in firstname seq. Skipped with opt_comparative as the
	       # length differences are already reflected in the relative score
#	       ( $opt_comparative  || $bit[5] >= $unalileft2-20 && $endUn >= $unaliright2-20 ) ) {
	    # is the first (query) also contained in the second (target)
	    if ( $begUn <= $unaligned && $endUn <= $unaligned
		&& $begUn + $endUn <= $totalunaligned ) { # &&
#		 ( $opt_comparative  || $bit[5] <= $unalileft2+20 && $endUn <= $unaliright2+20 ) ) {
	      # with opt_comparative the alignments lines that are read
	      # in here are ordered from high to low the script stops
	      # when $relscore gets below $opt_comparative and if the 
	      # opposite pair had a higher score, we would have looked at it already
	      if ($okaytoremove1 
		  && ($keep{$secname} 
		      # can only get here if !$keep{$firstone} ($goahead is only set to 1 if not both of them are on the keep list
		      || !$okaytoremove2 # same: only get here when $okaytoremove1
		      || $opt_table && ( $qualscore{$nopound2name} > $qualscore{$nopound1name} ||
					 $qualscore{$nopound2name} == $qualscore{$nopound1name} && $score2over1 )
		      || ( $opt_ambigfilter && 
			   (!$opt_table || $qualscore{$nopound2name} == $qualscore{$nopound1name}) &&
			   (!$ambigfrac{$secname} || $ambigfrac{$firstname} > $ambigfrac{$secname}) )
		      || !$opt_table && !$opt_ambigfilter && $score2over1 ) ) {
		&removeone($secname,$firstname,"1gone2 $_",$subfrag);
	      } elsif ($remove2isfine) {
		$entry{$secname} = "";
		&removeone($firstname,$secname,"2gone1 $_",$subfrag);		  
	      }
	    } elsif ($remove2isfine) {
	      # target contained in query, but not vice versa
	      &removeone($firstname,$secname,"2gone2 $_",$subfrag);
	    }
	  } elsif (!$keep{$firstname} && $okaytoremove1 &&
		   $begUn <= $unalileft2 + 20 && $endUn <= $unaliright2 + 20 &&
		   $begUn <= $unaligned && $endUn <= $unaligned &&
		   $begUn + $endUn <= $totalunaligned ) {
	    &removeone($secname,$firstname,"1gone2 $_",$subfrag);
	  } elsif ($opt_justhalf) {
	    if ($opt_table) {
	      # no dying, but will cause a (harmless ?) error a few lines down
	      warn "No qual for $firstname\n" if !defined $qualscore{$nopound1name};
	      warn "No qual for $secname\n" if !defined $qualscore{$nopound2name};
	    }
	    if (!$keep{$firstname} && $okaytoremove1 
		&& $match1frac > $match2frac && $match1frac > 0.5
		&& (!$opt_table || $qualscore{$nopound1name} < 2* $qualscore{$nopound2name})
		&& ($begUn <= $unalileft2 + 20 && $begUn <= $unaligned
		    || $endUn <= $unaliright2 + 20 && $endUn <= $unaligned) ) {
	      &removeone($secname,$firstname,"1goneH $_",$subfrag);
	    } elsif (!$keep{$secname} && $okaytoremove2
		     && $match2frac > $match1frac && $match2frac > 0.5
		     && (!$opt_table || $qualscore{$nopound2name} < 2*$qualscore{$nopound1name})
		     && ($begUn >= $unalileft2 - 20 && $unalileft2 <= $unaligned
			 || $endUn >= $unaliright2 - 20 && $unaliright2 <= $unaligned) ) {
	      &removeone($firstname,$secname,"2goneH $_",$subfrag);
	    }
	  } elsif ($merge && $bit[8] ne 'C') { # slight danger of not merging sequences simply named "C".
	    my $newcut = $div;
	    $newcut = 2 unless $div < 2 || $fragsoverlap;
	    if ($firstname eq $secname && $firstname ne $bit[8] && $firstname ne $bit[9]) {
## to avoid joining sequences head to tail (in a circle)
	      if ($endUn < 25 && $unalileft2 < 25) {
		($seq{$secname} = $entry{$secname}) =~ s/^>.*\n([A-Za-z]+)/$1/;
		print OUTINF "2mergedwith1 $_";
		$entry{$firstname} .= substr($seq{$secname},$bit[10]);
		$entry{$secname} = "";
	      } elsif ($unaliright2 < 25 && $begUn < 25) {
		($seq{$firstname} = $entry{$firstname}) =~ s/^>.*\n([A-Za-z]+)/$1/;
		print OUTINF "1mergedwith2 $_";
		$entry{$secname} .= substr($seq{$firstname},$bit[6]);
		$entry{$firstname} = "";
	      }
	    } elsif ($bit[1] < $newcut && $bit[2] < $newcut && $bit[3] < $newcut) {
	      if ($endUn < 2 && $unalileft2 < 3) {
		($seq{$secname} = $entry{$secname}) =~ s/^>.*\n([A-Za-z]+)/$1/;
		print OUTINF "2mergedwith1 $_";
		$entry{$firstname} .= substr($seq{$secname},$bit[10]);
		$entry{$secname} = "";
	      } elsif ($unaliright2 < 2 && $begUn < 2) {
		($seq{$firstname} = $entry{$firstname}) =~ s/^>.*\n([A-Za-z]+)/$1/;
		print OUTINF "1mergedwith2 $_";
		$entry{$secname} .= substr($seq{$firstname},$bit[6]);
		$entry{$firstname} = "";
	      }
	    }
	  }
	}
      } # if ($goahead)
    } # if cross_match line
  } # while reading crossmatchfile compared to self, masklevel 101 
} # unless $opt_nocrossmatch
close OUTINF;
open (OUT, ">$file.selection");
foreach my $key (sort { $orderN{$a} <=> $orderN{$b} } keys %entry) {
  if ( $exactas{$key} ) {
    # if the duplicates had the same classification, no reason to show that each time
    if ($key =~ /(\#\S+)/) {
      my $tag = $1;
      $exactas{$key} =~ s/$tag$//;
      $exactas{$key} =~ s/$tag,/,/g;
    }
    $exactas{$key} =~ s/,$//;
    $exactas{$key} =~ s/,,/,/g;
    $entry{$key} =~ s/\(a.k.a. \S+\)//; # remove previous exacta data (absorbed above)
    $entry{$key} =~ s/^(>.+?)\n/$1 \(a.k.a. $exactas{$key}\)\n/;
  }
  if ( $mention{$key} ) {
    if ($key =~ /(\#\S+)/) {
      my $tag = $1;
      $mention{$key} =~ s/$tag$//;
      $mention{$key} =~ s/$tag,/,/g;
    }
    $mention{$key} =~ s/,$//;
    $mention{$key} =~ s/,,/,/g;
    $entry{$key} =~ s/^(>.+?)\n/$1 Kept in favor of $mention{$key}\n/;
  }
  chomp $entry{$key}; # sometimes there is no \n at end
  print OUT "$entry{$key}\n" if $entry{$key};
}
close OUT;
unlink "tempxmatch.stderr" if -e "tempxmatch.stderr";
unlink "$query.log" if -e "$query.log";

sub getdata {
  # combines extra material on the fasta line of the two sequences being merged
  my ($data1,$data2) = @_;
  if ($data1 && !$data2) {
    return " $data1";
  } elsif (!$data1 && $data2) {
    return " $data2";
  } elsif ($data1 && $data2) {
    if ($data1 eq $data2 || $data1 =~ /$data2/) {
      # somehow $data1 =~ /$data2/ is false even if $data1 equals $data2
      return " $data1";
    } elsif ($data2 =~ /$data1/) {
      return " $data2";
    } else {
      return " $data1 $data2";
    }
  } else {
    return "";
  }
}

sub removeone {
  my ($stay,$go,$txt,$keepsubfrag) = @_;
  if ($keepsubfrag) {
    print OUTINF "$keepsubfrag";
  } else {
    warn "No name for \$go with $stay \n" unless $go;
    $mention{$stay} .= "$go,"; 
    $mention{$stay} .= "$mention{$go}," if $mention{$go};
    $entry{$go} = "";
    print OUTINF "$txt";
  }
}

sub checkifokay {

  # Check if it is okay to toss out a sequence if previously non-identical
  # sequences have been removed because of it (stored in the hash %mention. It
  # only checks if all the previously removed seqs have either a sufficiently
  # high relative score (with option -comparative) or has a divergence level
  # below the -d cutoff versus the other sequence. It is possible that too much
  # of a removed sequence would have been left unaligned to normally pass,
  # especially without the comparative option

  my $ok;
  my ($quer,$targ) = @_;
  $mention{$quer} =~ s/,,/,/g;
  my @tossees = split ',',$mention{$quer};
  foreach my $tossee (@tossees) {
    my $pair = "$tossee"."_$targ";
    if ($opt_comparative && !$relpairsc{$pair} ||
	!$opt_comparative && !$pairdiv{$pair}) {
        # The rejected sequence are listed on the fasta line of the "winner"
        # lack classification if it was the same as or was less specific than
        # that of the sequence it was replaced by (i.e. the query here).
      $tossee =~ /(\#\S*)$/;
      my $tclass = $1;
      $quer =~ /(\#\S+)$/;
      my $qclass = $1;
      if ($tclass && $qclass && $tclass ne $qclass) {
	$tossee =~ s/$tclass/$qclass/;
      } elsif (!$tclass && $qclass) {
	$tossee .= $qclass;
      }
      $pair = "$tossee"."_$targ";
    }
    if ($opt_comparative) {
      if (!$relpairsc{$pair}) {
	# I have not tossed any alignments but there are a number of reasons this may happen:
	# - the score of tossee vs targ was below cutoff (and just above cutoff vs $quer)
	# More commonly:
	# - $tossee's match to $quer doesn't fully overlap the $targ vs $quer
	# match (meaning $targ is smaller than $quer). These alignments are
	# doomed anyway as the comparative score of $quer vs $targ should be too
	# small. I should see if I can stop this earlier
	# - The problematic one:
	# The file that one tries to deredundate was created by an earlier round
	# of deredundification. The fasta lines contain information
	# "Kept in favor of DFxx" that is picked up by this script and the entry DFxx
	# is not in the fasta file so not compared to the target in this run.
	# One solution is to specifically add from the previous round of
	# comparisons from the .vse file those alignment-lines (no full
	# alignments necessary) that it complains about to the current .vse file
	# and redo the analysis using -r.
	# This becomes impossible when the sequences were in different clusters 
	# (if one applied a quick clustering). For now, the solution would be to
	# grab those  from the original (completely unduplicated) fasta files and
	# do the match there and then... That's a pain.
	$ok = 0; 	
	print STDERR "No \$relpairsc for $tossee (tossee of $quer) vs $targ\n";
      } elsif ( $relpairsc{$pair} >= $opt_comparative) {
	$ok = 1;
      } else {
	$ok = 0;
	last;
      }
    } else {
      if (defined $pairdiv{$pair} && $pairdiv{$pair} <= $div) {
	$ok = 1;
      } else {
	$ok = 0;
	last;
      }
    }
  }
  return $ok;
}

sub onpath {
  my $prog = shift;
  foreach my $dir (split /:/, ($ENV{'PATH'} || "")) {
    next unless length $dir;
    return "$dir/$prog" if -x "$dir/$prog" && ! -d "$dir/$prog";
  }
  return;
}

sub needonpath {
  my @missing = grep { !&onpath($_) } @_;
  return unless @missing;
  die "Cannot find " . join(", ", @missing) . " on your PATH.\n" .
      "cross_match comes with phrap; rmblast.pl with RepeatModeler; " .
      "rmblastn and makeblastdb with rmblast.\n";
}
