#!/usr/bin/perl -w

# A shortcut to run cross_match with 
# - shorter options to fill in
# - control of noise in STDERR and STDOUT 
# - recalculates mismatch level as most people would expect to see it
# - option to display score per base relative to maximum (in %)
# - option to discount one transition per CpG site (only implemented with above option)

#20200415: addition of -p option and some straightening up of the output
#20210401: added  -cg option; fixed max score expectations for sequences with ambiguous bases in it
#20210422: -x option now delivers a masked fasta file with sequences on one line, rather than in blocks of 60 bp
#20210527: some more quotemetas to deal with targets that have pipes in them
#20211203: Now reads matrix; fixes and improvements in -p score calculation. 
#20220511: Included rmblast.pl option (-blast), -forward option. Cleaned up output even more.
#20220517: Allow multiple subjct files (temporarily combined in one file)
#20220520: Fixed some behavior with -zip option; allowed adjustment of -score with -zip
#20220921: Sped up CpG adjustments to score
#20230915: fixed /u1/home -> /home directories and >& redirects
#20260902: matrices moved to Matrices/ next to the script; cross_match and
#          rmblast.pl now have to be on the PATH; fixed -zip deleting the
#          cross_match output it was handed

# TO DO
# - implement -cg outside -p
# - implement -forward to create easier to parse tables

# One feature of cross_match that is neither copied by Xmatch.pl nor
# rmblast.pl is that cross_match can take multiple subject files:
# "cross_match file1 file2 file3"
# results in a comparison of file1 to file2 and file3 combined


use strict;
use Getopt::Long;
use FindBin;

# Matrices ship with this script, under Matrices/ in the installation directory.
my $matrixdir = "$FindBin::RealBin/Matrices";


my $USAGE = "Usage: Xmatch.pl <fasta seq1> (<fasta seq2>)

option           equivalent cross_match option 
-a(lignments)    -alignments
-ba(ndwidth)     -bandwidth (default 14)
-d(el_gap_ext)   -del_gap_ext; with this option -ins_gap_ext is set to (del_gap_ext - 1) for a higher penalty
-e(xtension)     -gap_ext ;  sets both -del_gap_ext and -ins_gap_ext (default -5)
-g(ap_init)      -gap_init (default -25)
-l(evel)         -masklevel (default 80)
-ma(trix)        -M / -matrix (default Matrices/crossmatch/ctools.matrix, or
                 Matrices/ncbi/nt/ctools.matrix with -blast)
                 if a number 14,18,20 or 25 is given ##p41g.matrix is chosen
                 if given 14p35 (37,39,41,43,45,47,49,51,53) #####g.matrix is chosen
                 otherwise give a path, or a file name in the matrix directory
-mi(nmatch)      -minmatch (default 14)
-r(aw)           -raw (default off, i.e. score is complexity adjusted) 
-s(core)         -minscore (default 200)
-w(ord_raw)      -word_raw (default off, i.e. length of seeds is higher with lower complexity)
-x               -screen (default off)

-bl(ast)         use rmblast instead of cross_match (options are the same) 
-o(riginal)      keep the original cross_match mismatch level (mismatches/length_of_query) 
                 By default Xmatch.pl recalculates the mismatch level as mismatches/aligned_bases
-q(uiet)         STDERR of cross_match gets redirected to a file \"tempxmatch.stderr\" 
-n(ab)           Only prints the summary lines and adds a \"+\" for forward match to give constant nr of columns
-zip <xmoutput>  Skip the cross_match part and take the indicated cross_match output or Xmatch.pl .rawxmout output
                 (e.g. to increase the minimum score displayed (-score) or to add or skip -perbasescore, -cg, or -nab)

   # Not reliable with -blast option, until rmblast.pl has a -raw option: 
-p(erbasescore)  Adds a first column reporting the raw score per base of the query 
                 NOTE1: this forces the -raw option, so only use this on established matches, and the -cg option
                 NOTE2:  unless a -matrix option is used, this mode will use symmetricNscores0matrix 
   #  For now only works within a -p run.
-cg              Only works and effects a -p run right for now. Subtracts from mismatch level and adjusts up scores for
                 one transition/ambiguous mismatch per CG in either consensus and treats TG/CA pairs as one transition.
                 Forces -a (as alignments are needed for the adjustment)
";

my $commandline = join " ", $0, @ARGV;
print "$commandline\n"; #better print to STDERR or to some little .log file? cross_match prints its to STDOUT

my @opts = qw (alignments bandwidth=s blast cg del_gap_ext=i extension=i gap_init=i level=i matrix=s minmatch=i nab original perbasescore quiet score=i raw word_raw x zip=s);
our ($opt_alignments, $opt_bandwidth, $opt_blast,$opt_cg, $opt_del_gap_ext, $opt_extension, $opt_gap_init, $opt_level, $opt_matrix, $opt_minmatch, $opt_nab, $opt_original, $opt_perbasescore, $opt_quiet, $opt_score, $opt_raw, $opt_word_raw, $opt_x,$opt_zip);

unless (&GetOptions(@opts) && $ARGV[0]) {
  warn "\n$USAGE";
  exit (1);
}

die "$USAGE\n$ARGV[0] is not a readable file\n" unless -s $ARGV[0] || $opt_zip;
my $fastafile1 = shift if $ARGV[0];
my $fastafile2 = "";
my $selfcomparison;
my $subfilecnt;
my $tempsubjfile;
while ($ARGV[0]) {
  if (-s $ARGV[0]) {
    if ($subfilecnt) {
      my $extrafile = shift;
      my $efname = $extrafile;
      $efname =~ s/^\S+\///; # remove any directory structure
      if ($tempsubjfile) {
	$tempsubjfile .= "_$efname";         
      } else {
	$tempsubjfile = $fastafile2;
	$tempsubjfile =~ s/^\S+\///; # remove any directory structure
	$tempsubjfile = "tempcombinedfiles_$tempsubjfile\_$efname";
      }
      system "cat $fastafile2 $extrafile > $tempsubjfile";
      $fastafile2 = $tempsubjfile;
    } else {
      $fastafile2 = shift;
    }
    ++$subfilecnt;
  } else {
    die "$USAGE\n$ARGV[0] is not a readable file\n"
  }
}
if ($opt_blast && !$fastafile2) {
  # when cross_match notices that there is not subject file, it assumes
  # masklevel 101 (it literally sets masklevel to 101) and does not
  # show self matches (which it does show in "cross_match file file")
  # Logical as there is no use for an output of self matches only.
  # rmblast does not work without a subject file. Since this is a
  # common use, cross_match's behavior is copied by:
  $fastafile2 = $fastafile1;
  $opt_level = 101;
  $selfcomparison = 1;
  # currently the outputs are different in that with rmblast you get the self matches
  # while cross_match does not report these
}

# cross_match and rmblast.pl have to be on the PATH. rmblast.pl comes with
# RepeatModeler and calls rmblastn and makeblastdb.
unless ($opt_zip) {
  if ($opt_blast) {
    &needonpath("rmblast.pl","rmblastn","makeblastdb");
  } else {
    &needonpath("cross_match");
  }
}

my $matrix;
# default matrices. R vs A/G/R and Y vs C/T/Y do not score symmetrically
# in ctools.matrix.
if ($opt_blast) {
  $matrix = "$matrixdir/ncbi/nt/ctools.matrix";
} else {
  $matrix = "$matrixdir/crossmatch/ctools.matrix";
}
my $matdir = $matrix;
$matdir =~ s/ctools.matrix$//;

my $command; # will contain the cross_match or rmblastn commandline
if ($opt_zip) {
  die "$USAGE\n$opt_zip is not a readable file\n" unless -s $opt_zip;
  my $sawxmline;
  open (TEST, $opt_zip);
  while (<TEST>) {
    if (/^cross_match/ && !$sawxmline) {
      unless (/alignments/) {
	die "With -cg or -perbasescore and -zip you need a cross_match output with alignments as input.\n" if $opt_cg || $opt_perbasescore;
      }
      $command = $_;
      if ($opt_score) {
	$command =~ s/-score \d+/-score $opt_score/ || $command =~ s/\n$/-score $opt_score\n/;
      }
      # extract the matrix used; other parameters do not affect recalculations when changing -p, -cg or -f options
      if (/-matrix\s+(\S+)/) {
	$matrix = $1;
      } else {
	$matrix = "$matrixdir/crossmatch/defaultcrossmatchmatrix";
      }
      ++$sawxmline;
    } elsif ( /^Query file\(s\):\s+(\S+)/) {
      # not sure yet if the query and subject files are needed for anything 
      $fastafile1 = $1;
    } elsif (/^Subject file\(s\):\s+(\S+)/) { # not always there
      $fastafile2 = $1;
      last;
    } elsif (/^Presumed sequence/) {
      last;
    } elsif ( $. > 10) {
      last;
    }
  }
  die "With -zip use an unmodified cross_match output or the .rawxmout output of Xmatch.pl\n" unless $sawxmline;
  close TEST;
}

if ( $opt_perbasescore ) {
  $opt_cg = 1;
} elsif ($opt_cg) {
  die "The option -cg has only implemented within the option -perbasescore for now.\n";
}

my $alignments = "";
$alignments = "-alignments" if $opt_alignments || $opt_cg;

my $bandwidth = ""; 
if (defined $opt_bandwidth) {
  if ($opt_bandwidth < 0) {
    die "$USAGE\nbandwidth (-b) needs to be positive or 0 for ungapped alignments\n"; 
  } else {
    $bandwidth = "-bandwidth $opt_bandwidth";
  }
}
my $gapext = "-gap_ext -5";
if (defined $opt_del_gap_ext) {
  die "$USAGE\nEither choose gap_ext (-e) or del_gap_ext (-d), not both" if $opt_extension;
  $opt_del_gap_ext = -$opt_del_gap_ext if $opt_del_gap_ext > 0;
  if ($opt_del_gap_ext == 0 ) {
    die "$USAGE\ndel_gap_ext (-d) should be negative\n";
  } else {
    my $ins = $opt_del_gap_ext - 1;
    $gapext = "-del_gap_ext $opt_del_gap_ext -ins_gap_ext $ins";
  }
}
if (defined $opt_extension) {
  $opt_extension = -$opt_extension if $opt_extension > 0;
  if ($opt_extension == 0) {
    die "$USAGE\ngap_ext (-e) should be negative";
  } else {
    $gapext = "-gap_ext $opt_extension";
  }
}
my $gapinit = "-gap_init -25";
if (defined $opt_gap_init) {
  $opt_gap_init = -$opt_gap_init if $opt_gap_init > 0;
  if ($opt_gap_init == 0 ) {
    die "$USAGE\ngap_init (-g) should be negative \n";
  } else {
    $gapinit = "-gap_init $opt_gap_init";
  }
}

my $masklevel = ""; # cross_match defaults to -masklevel 80
if (defined $opt_level) {
  if ($opt_level < 0  && $opt_level > 101) {
    die "$USAGE\nmasklevel (-l) is outside the range of 0-101";
  } else {
    $masklevel = "-masklevel $opt_level";
  }
}

unless ($opt_zip && $matrix) {
  if ($opt_matrix) {
    if ($opt_matrix =~ /^\d+$/) {
      if ($opt_matrix == 14 || $opt_matrix == 18 || $opt_matrix == 20 || $opt_matrix == 25) {
	$matrix = "$matdir$opt_matrix"."p41g.matrix";
      } else {
	die "Can only choose from 14, 18, 20 or 25 for matrix divergence levels\n";
      }
    } elsif ($opt_matrix =~ /^\d{2}p\d{2}g?$/) {
      $matrix = "$matdir$opt_matrix"."g.matrix";
      $matrix =~ s/gg\.matrix/g\.matrix/;
      die "$matrix does not seem to exist at $matdir\n" unless -s $matrix;
    } elsif (-s $opt_matrix) {
      $matrix = "$opt_matrix";
      if ($opt_blast) {
	$matrix =~ s/\S+\//$matdir/;
      }
    } else {
      $matrix = "$matdir$opt_matrix";
      die "Cannot find $opt_matrix in current directory or at $matdir\n" unless -s $matrix;
    }
  } elsif ($opt_perbasescore) {
    # $matdir already points at the rmblast matrices when -blast is set
    $matrix = "$matdir"."symmetricNscores0matrix";
  }
}

# read in matrix
my @base;
my %matscore;
my $j = 0;
open (INM, "$matrix") || die "File $matrix does not exist\n";
while (<INM>) {
  next if /^\#?\s*FREQS/;
  if (/^\s*[A-Za-z]\s+[A-Za-z]/) {
    tr/a-z/A-Z/;
    @base = split;
  } elsif (/^([A-Za-z])?\s*\-?\d/) {
    my @score = split;
    shift @score if $score[0] =~ /[A-Za-z]/;
    for (my $i = 0; $i <= $#base; ++$i) {
      # tested with 14p41g.matrix : A in top sequence $i (genomic
      # copy) aligned to G in bottom sequence $j (e.g. consensus) gets
      # a penalty of -7 , vice versa -11
      $matscore{$base[$i]}{$base[$j]} = $score[$i];
    }
    ++$j;
  }
}

if ($opt_perbasescore) {
  # establish score of base to self if not present in matrix (cross_match treats them like N)
  foreach my $b ('A'..'Z') {
    $matscore{$b}{$b} = $matscore{'N'}{'N'} unless defined $matscore{$b}{$b};
  }
  # test for symmetry 
  my $wtext = "\nThis may lead to unexpected output with the -perbasescore option.
You could symmetrise your matrix or use the default matrix.\n";
  for (my $i = 0; $i <= $#base; ++$i) {
    for (my $j = 0; $j <= $#base; ++$j) {
      if ( $matscore{$base[$i]}{$base[$j]} != $matscore{$base[$j]}{$base[$i]} ) {
	die "$matrix is not completely symmetric (seq1 vs seq2 could score
differently then seq2 vs seq1). For example:
$base[$i] to $base[$j] scores $matscore{$base[$i]}{$base[$j]} and $base[$j] vs $base[$i] scores $matscore{$base[$j]}{$base[$i]}$wtext";
      }
    }
  } 
  if ( $matscore{'G'}{'G'} != $matscore{'C'}{'C'} ) {
    die "$matrix appears strand specific. 
It scores G:G different ($matscore{'G'}{'G'}) than C:C ($matscore{'C'}{'C'})$wtext";
  } elsif ( $matscore{'A'}{'A'} != $matscore{'T'}{'T'} ) {
    die "$matrix appears strand specific.
It scores A:A different ($matscore{'A'}{'A'}) than T:T ($matscore{'T'}{'T'})$wtext";
  } 
}
$matrix = "-matrix $matrix"; 
my $CtoTpenalty = $matscore{'C'}{'C'} - $matscore{'C'}{'T'}; # used in calccgadjust sub over and over

my $minmatch = "";
if ($opt_minmatch) {
  if ($opt_minmatch < 1) {
    die "$USAGE\nminmatch (-m) needs to be a positive number";
  } else {
    $minmatch = "-minmatch $opt_minmatch";
  }
}

my $minscore = "-minscore 200";
if ($opt_score) {
  if ($opt_score < 1) {
    die "$USAGE\nminscore (-s) needs to be a positive number\n";
  } else {
    $minscore = "-minscore $opt_score";
  }
}

# skippable with $opt_zip, but need to rearrange order;
my $quiet = "";
$quiet = 1 if $opt_quiet;
my $raw = "";
my $wordraw= "";
my $screen = "";
if ($opt_blast) {
  print "\nWARNING: The options -raw, -wordraw or -screen can not yet be used with the rmblastn engine (-blast option).
The -raw, -wordraw and -screen options are ignored in this run.\n\n" if $opt_raw || $opt_word_raw || $opt_x;
  print "\nWARNING: Lacking a -raw option, the relative score produced by the -perbasescore option with the rmblastn engine is dependent on the complexity of the sequence. The relative score of low complexity sequences will be lower then they should be.\n\n" if $opt_perbasescore;
} else {
  $raw = "-raw" if $opt_raw || $opt_perbasescore;
  $wordraw = "-word_raw" if $opt_word_raw;
  $screen = "-screen" if $opt_x;
}

# also skippable with $opt_zip, except for the %ambigadj data; perhaps write that to the rawxmout file?

# check for entries with identical names
# count ambiguous bases in query to calculate proper maximum score for $opt_perbasescore
my %ambigadj;
my $filecounter;
my $foundadupname;
foreach my $bestand ($fastafile1, $fastafile2) {
  # "bestand" is a Dutch word for file
  last unless $bestand; # i.e. when $fastafile2 does not exist
  if ($selfcomparison) {
    $fastafile2 .= '.modified' if $foundadupname; # $fastafile1 name has been changed to "$fastafile1.modified"
    last;
  }
  my %seenthat;
  my ($seqname,$seqstr,$allwhatsinbestand);
  my $modfile = "$bestand.modified";
  open (IN, "$bestand");
  while (<IN>) {
    if ( /^>/ ) {
      if ($opt_perbasescore &&
	  !$filecounter && # adjustments are only made for the first file
	  $seqstr ) {
	&calcambigadj($seqname,$seqstr);
	$seqstr = "";
      }
      $_ =~ /^>(\S+)/;
      $seqname = $1;
      if ($seenthat{$seqname}) {
	my $newseqname = "$seqname.\[$seenthat{$seqname}\]";
	my $occ = $seenthat{$seqname} + 1;
	warn "$seqname occurs $occ times in $bestand. Nr$occ renamed to $newseqname in output and $modfile\n";
	s/^>$seqname/>$newseqname/;
	$seqname = $newseqname; # to be fed to calcambigadj
	$foundadupname = 1;
      }
      ++$seenthat{$seqname};
    } elsif (/\w/) {
      $seqstr .= $_;
    }
    $allwhatsinbestand .= $_;
  }
  &calcambigadj($seqname,$seqstr) if $opt_perbasescore && $seqstr;
  if ($foundadupname) {
    open (OUTMOD, ">$modfile");
    print OUTMOD "$allwhatsinbestand";
    close OUTMOD;
    if ($bestand eq $fastafile1) {
      $fastafile1 = $modfile;
    } elsif ($bestand eq $fastafile2) {
      $fastafile2 = $modfile;
    }
  } 
  $allwhatsinbestand = "";
  ++$filecounter;
}

unless ($opt_zip) {
  # check if the subject file blast database is up to date
  if ($opt_blast) {
    my $ndbfile = "$fastafile2.ndb";
    if ( -s $ndbfile && (stat($ndbfile))[9] < (stat($fastafile2))[9] ) {
      unlink ($ndbfile, "$fastafile2.nhr", "$fastafile2.nin", "$fastafile2.nog", "$fastafile2.nos", "$fastafile2.not", "$fastafile2.nsq", "$fastafile2.ntf", "$fastafile2.nto");
    }
  }

  # execute cross_match command
  $command = "$fastafile1 $fastafile2 $alignments $bandwidth $gapext $gapinit $masklevel $matrix $minmatch $minscore $raw $wordraw $screen";
  if ($opt_blast) {
    $command = "rmblast.pl $command";
  } else {
    $command = "cross_match $command";
  }
}

print "$command\n" if $opt_nab;# otherwise already printed to STDOUT by both rmblast.pl and cross_match
my $rawoutfile = "$fastafile1"."_vsit.rawxmout";
my $short2 = $fastafile2;
$short2 =~ s/^\S+\///;
$rawoutfile = "$fastafile1"."_vs$short2.rawxmout" if $fastafile2;
if ($opt_zip) {
  $rawoutfile = $opt_zip;
} elsif ($quiet) {
  system "$command 1> $rawoutfile 2> tempxmatch.stderr";
} else {
  system "$command > $rawoutfile";
}

# print sequences in .screen file on one line each
if ($screen) {
  my $xdfile = "$fastafile1.screen";
  if ($opt_zip) {
    warn "Currently, no (new) masked file $xdfile is created when skipping the cross_match / rmblastn step with the -zip option\n\n";
  } else {
    # cross_match prints out a.screen file in 50 letter blocks
    # something to do with limited memory of computers in the past?
    # I much prefer sequences to be on one line.
    my $random = rand(999999);
    my $tempxdfile = "extraordinarilytemporaryfile$random";
    open (OUT, ">$tempxdfile");
    open (IN, "$xdfile");
    my $linecnt = 0;
    while (<IN>) {
      if (/^>/) {
	if ($linecnt) {
	  print OUT "\n$_";
	} else {
	  print OUT;
	}
      } else {
	chomp;
	print OUT;
      }
      ++$linecnt;
    }
    print OUT "\n";
    close OUT;
    rename $tempxdfile,$xdfile;
  }
}

# clean up the temporary file that is a concatenation of subject files
if ($tempsubjfile) {
  system "rm $tempsubjfile*";
}


# I suppose this is faster than twice reading a file in line by line
open (IN, "$rawoutfile");
my @lines = <IN>;
close IN;

my $print = 1;

if ($opt_perbasescore) {
  print "The first column is the (raw) score divided by the (raw) score of the query to
itself times 100. If a non-symmetrical matrix was optioned for the alignment,
this fraction will depend on the query sequence composition.\n\n"
}

# run through lines to make adjustments with $opt_cg
my (%addtoscore,%nrmismatches) = ();
if ($opt_cg) {
  warn "Making adjustments for CpG sites\n" unless $opt_quiet;
  my ($xmline,$seq1,$mutline,$seq2,$turn,$leftspace) = ();
  foreach ( @lines ) {
    if ( /^\s*\d\d.+\(\d+\)/ ) {
      if ($xmline) {
	($addtoscore{$xmline},$nrmismatches{$xmline}) = &calccgadjust($seq1,$mutline,$seq2);
	$seq1 = "";
	$mutline = "";
	$seq2 = "";
	$turn = 0; #should already be anyway
      }
      $xmline = $_;
    } elsif (/^C?\s+\S+\s+\d+\s([A-Z\-]+)\s\d+/) {
#e.g.
#   AciJub-1.56#LIN         1 ATGGAATATTACTCGGC-ATCAAAAAGAATGAAATCTTGCCATTTGCAAC 49
# or
# C CrasThon-1.62#D       193 TAGGCTGACCAGACGTACCATTTTAAATGGGACAGAAAACAACGGTACTG 144
      if ($turn) {
	$seq2 .= $1;
	$turn = 0;
      } else {
	$seq1 .= $1;
	$turn = 1;
      }
      unless ($leftspace) {
	/^(C?\s*\S+\s+\d+\s+)/;
	$leftspace = length $1;
      }
    } elsif ($turn) {
      $mutline .= substr($_,$leftspace);
      chomp $mutline;
    }
  }
  if ( $xmline ) {
    ($addtoscore{$xmline},$nrmismatches{$xmline}) = &calccgadjust($seq1,$mutline,$seq2);
  }
}

my (%max0,%max5,%max6,%max7,%max9,%max10,%max11,%max12);
# to hold the max length of score, leftinquery, subject_name, and leftinsubj
my @newlines;
my ($skipselfalignment,$skipforlowscore);
 warn "Adjusting format\n" unless $opt_quiet;
foreach ( @lines ) { 
  # the following eliminates the binary-like stuff that appears
  # when cross_match encounters a base in the matrix it doesn't 
  # know about (for me usually "Z")
  if ( /^Num. pairs/ || /^Discrepancy summary/) {
    $print = 0;
  } elsif ( !$print && /\(\d+\)/ ) {
    $print = 1;
  }
  if ($print) {
    if ( /\(\d+\)/ ) {
      my $cmline = $_;
      my @bit = split;
      next unless $bit[11];
      if ($selfcomparison && $bit[4] eq $bit[8] 
	  && $bit[5] eq $bit[9] && $bit[7] eq $bit[11]) {
	# rmblast.pl run of file against itself (no subject file indicated)
	# cross_match does not report selfmatches, so let's have rmblast.pl also not do this
	$skipselfalignment = 1;
      } else {
	$skipselfalignment = 0;
	if ( $addtoscore{$_} ) { # only happens wit $opt_cg
	  # also adjusting the mismatch level in column2/bit[1] by
	  # considering CpG to TG or CA not as a mismatch each one of
	  # these represents 10 points in the $addtoscore
	  $bit[0] += $addtoscore{$_};
	  s/^(\s*)\d+/$1$bit[0]/;
	}
	if ($opt_zip && $opt_score) {
	  if ($opt_score > $bit[0]) {
	    $skipforlowscore = 1;
	  } else {
	    $skipforlowscore = 0;
	  }
	}
	unless ($opt_original || $bit[1] eq '0.00') {
	  my $ql = $bit[6] - $bit[5] + 1;
	  my $mismatches = int ($ql * $bit[1] / 100 + 0.5); 
	  # proper rounding; if the number was rounded properly to start with
	  $mismatches = $nrmismatches{$cmline} if defined($nrmismatches{$cmline}); 
	  # only happens with $opt_cg; can be zero now
	  my $alignedbases = int ($ql * (1-$bit[3]/100) + 0.5);
	  my $realmmlevel = sprintf "%4.2f", (100*$mismatches/$alignedbases);
	  if ($realmmlevel < 10) {
	    s/^(\s*\d+)\s+$bit[1]/$1  $realmmlevel/;
	  } else {
	    s/^(\s*\d+)\s+$bit[1]/$1 $realmmlevel/;
	  }
	}
	# establish maximum width fields for each query name
	# empty lines separate different queries so it does
	# not look chaotic to have different widths for each
        $max0{$bit[4]} = length $bit[0] if !$max0{$bit[4]} || (length $bit[0]) > $max0{$bit[4]};
        $max5{$bit[4]} = length $bit[5] if !$max5{$bit[4]} || (length $bit[5]) > $max5{$bit[4]};
        $max6{$bit[4]} = length $bit[6] if !$max6{$bit[4]} || (length $bit[6]) > $max6{$bit[4]};
        $max7{$bit[4]} = length $bit[7] if !$max7{$bit[4]} || (length $bit[7]) > $max7{$bit[4]};
	splice @bit,8,0,'+' unless $bit[8] eq 'C' && $bit[12];
	$max9{$bit[4]} = length $bit[9] if !$max9{$bit[4]} || (length $bit[9]) > $max9{$bit[4]};
	$max10{$bit[4]} = length $bit[10] if !$max10{$bit[4]} || (length $bit[10]) > $max10{$bit[4]};
	$max11{$bit[4]} = length $bit[11] if !$max11{$bit[4]} || (length $bit[11]) > $max11{$bit[4]};
	$max12{$bit[4]} = length $bit[12] if !$max12{$bit[4]} || (length $bit[12]) > $max12{$bit[4]};
      }
      push @newlines, $_ unless $skipselfalignment || $skipforlowscore;
    } else {
      push @newlines, $_ unless $opt_nab ||$skipselfalignment || $skipforlowscore;
    }
  }
}
@lines = ();
my $lastq;
foreach ( @newlines ) {
  if ( /\(\d+\)/ ) {
    # Format space for clearer output. cross_match only adjusts the
    # distance between query and start position for the longest
    # length in the column
    my @bit = split;
    print "\n" if $opt_nab && (!$lastq || $lastq ne $bit[4]);
    my ($sp0,$sp5,$sp6,$sp7,$sp9,$sp10,$sp11,$sp12);
    my $sp1 = "";
    my $sp2 = "";
    my $sp3 = "";
    $sp0 = " " x ($max0{$bit[4]} - (length $bit[0]) );
    $sp1 = " " if $bit[1] < 10;
    $sp2 = " " if $bit[2] < 10;
    $sp3 = " " if $bit[3] < 10;
    $sp5 = " " x ($max5{$bit[4]}-(length $bit[5]));
    $sp6 = " " x ($max6{$bit[4]}-(length $bit[6]));
    $sp7 = " " x ($max7{$bit[4]}-(length $bit[7]));
    splice @bit,8,0,'+' unless $bit[8] eq 'C' && $bit[12];
    $sp9 = " " x ($max9{$bit[4]}-(length $bit[9]));
    $sp10 = " " x ($max10{$bit[4]}-(length $bit[10]));
    $sp11 = " " x ($max11{$bit[4]}-(length $bit[11]));
    $sp12 = " " x ($max12{$bit[4]}-(length $bit[12]));
    $bit[8] = " " if $bit[8] eq '+' && !$opt_nab;
    $_ = "$sp0$bit[0] $sp1$bit[1] $sp2$bit[2] $sp3$bit[3]  $bit[4]   $sp5$bit[5] $sp6$bit[6] $sp7$bit[7]  $bit[8] $bit[9]$sp9  $sp10$bit[10] $sp11$bit[11] $sp12$bit[12]\n";

    # adding the extra column with per base similarity score
    if ($opt_perbasescore) {
      $bit[7] =~ tr/()//d;
      # non-complexity adjusted score of the full sequence against itself
      my $maxsc = ($bit[6] + $bit[7])*$matscore{'A'}{'A'};
      # modify expected score vs itself given ambiguous bases in the sequence
      $maxsc -= $ambigadj{$bit[4]} if $ambigadj{$bit[4]};
      my $perbase = sprintf("%4.1f", 100*$bit[0]/$maxsc);
      if ($perbase >= 100) {
	print "100.  $_";
      } else {
	print "$perbase  $_";
      }
    } else {
      print "$_";
    }
    $lastq = $bit[4];
  } else {
    print "$_";
  }
}

unlink $rawoutfile unless $opt_zip; # with -zip this is the caller's own file

sub calccgadjust {
  my $sq1 = shift; # full aligned sequence (in upper case, with dashes and all) of query
  my $mutationline = shift; # with i v - and ?
  my $sq2 = shift; # full aligned sequence of target
  my $adjustment = 0; # has to set to zero

  my $mismatches = ( $mutationline =~ tr/iv// );
  # Note that this does not count question marks and most aligned
  # ambiguous bases are not counted as mismatches. This is what I want
  # in $opt_perbase , but may not be always right when calccgadjust is
  # used outside -perbase For example, K:A, W:S, M:G etc should always
  # be counted as a mismatch. (R:Y already gets a "v" by
  # cross_match). Those are some rare situations, I suppose.

  # are there CpGs in the aligned part of the query?
  my $sq1cp = $sq1;
  my $mutlinecp = $mutationline;
  my ($len1,$len2);
  while ($sq1cp =~ s/^(\S*?)C(\-*)G//) { # query has a CG pair; change it to cg
    $len1 = length $1; #distance to C of first remaining CG
    $len2 = $len1 + (length $2) + 1; # distance to G of that CG (incl the C)
    # (length $2) representing the number of dashes between C and G

    if ( $mutlinecp =~ /^.{$len1}i/ || $mutlinecp =~ /^.{$len2}i/ ) {  
      # cross_match indicated a transition; there will be only one adjustment for "ii" (i.e. CG->AT) 
      $adjustment += 4*$CtoTpenalty/5;
      --$mismatches;
      # $CtoTpenalty equals $matscore{'C'}{'C'} - $matscore{'C'}{'T'}
      # what this does is reduce the original difference of a match
      # and a mismatch to the total score to 1/5th of that difference

    } elsif ( $mutlinecp =~ /^.{$len1}\?\s/ || $mutlinecp =~ /^.{$len2}\?/ ) {
      # this is a CpG lined up to NG or CN in the second sequence;
      # since CpGs have a lot of variability NG or CN are common
      # outcomes The adjustment is small, because in the -perbase mode
      # we do not want to make CN and CG calls equally good
      $adjustment += $CtoTpenalty/3;
      # no mismatch subtraction as the question marks were not counted
    # (could do different adjustments when CpG options works outside perbase option)
    }
    $mutlinecp =~ s/^.{$len2}.//; #the extra period for the removed G
  }
  # are there CpGs in the aligned part of the target? 
  # doesn't matter if the query has a CpG there too, as no substitutions will be reported
  $mutlinecp = $mutationline;
  while ($sq2 =~ s/^(\S*?)C(\-*)G//) {
    my $len1 = length $1;
    my $len2 = $len1 + (length $2) + 1;
    if ( $mutlinecp =~ /^.{$len1}i/ || $mutlinecp =~ /^.{$len2}i/ )  {
      --$mismatches;
      $adjustment += 9*$CtoTpenalty/10;
      # difference with match is now 1/10 instead of 1/5 above, since we favor 
      # models with CpG calls over ones with CA and TG calls
      # no adjustment for NG or CN opposite CG in target, as Ns were
      # already adjusted (+9 or +10) when counting the ambiguous bases in the query
    }
    $mutlinecp =~ s/^.{$len2}.//;
  }
  # adjust for CA opposite TG
  # to speed things up, we skip looking for situations like:
#  yoh            1 GCATTGCATGGTTTAGATAGCTTCGCCATT---GTGTAGTGGTA 41
#                                                i---i          
#  yah            1 GCATTGCATGGTTTAGATAGCTTCGCCATCTTAATGTAGTGGTA 44
  # also as it seems less likely that this involved a CpG transition 
  # It could  if the insert in "yah" starts with an A or ends with a C,
  # but it appears that cross_match aligns the transitions on one side
  # of a gap of more than 1 bp:
#  yoh            1 GCATTGCATGGTTTAGATAGCTTCGCCAT---TGTGTAGTGGTA 41
#                                                ---ii          
#  yah            1 GCATTGCATGGTTTAGATAGCTTCGCCATCTTCATGTAGTGGTA 44
  # and
#  yah            1 GCATTGCATGGTTTAGATAGCTTCGCCATCACAATGTAGTGGTA 44
#                                                ii---          
#  yoh            1 GCATTGCATGGTTTAGATAGCTTCGCCATTG---TGTAGTGGTA 41
  # but:
#  yah            1 GCATTGCATGGTTTAGATAGCTTCGCCATCAATGTAGTGGTA 42
#                                                i-i          
#  yo             1 GCATTGCATGGTTTAGATAGCTTCGCCATT-GTGTAGTGGTA 41

# For CG vs CA/TG 
#  AciJub-6.1674#L         1 GATGGCGGAACAGCATGGAAGTTTTTTGC--GTCTCTCGTCCATGAAATA 48
#                                                       vi--     v             
#  ParHer-4.26#LIN        12 GATGGCGGAACAGCATGGAAGTTTTTTTTCTGTCTCGCGTCCATGAAATA 61
# so I'll keep the check for multi-base gaps above, 

  $sq1cp = $sq1;
  $mutlinecp = $mutationline;
  while ($sq1cp =~ s/^(\S*?)C(\-?)A//) { # *? for the shortest distance
    my $len = length $1;
    my $gap = $2;
    if (!$gap && $mutlinecp =~ /^.{$len}ii/ 
	 || $gap && $mutlinecp =~ /^.{$len}i\-i/ ) {
      --$mismatches;
      $adjustment += 9*$CtoTpenalty/10;
      # same adjustment as for CG opposite CA or TG, since here one transition
      # represents a CpG change, the other a normal transition
    }
    $len += 2; # for the C and A
    ++$len if $gap;
    $mutlinecp =~ s/^.{$len}//;
  }
  # adjust for TG opposite CA

  $sq1cp = $sq1;
  $mutlinecp = $mutationline;
  while ($sq1cp =~ s/^(\S*?)T(\-?)G//) {
    my $len = length $1;
    my $gap = $2;
    if (!$gap && $mutlinecp =~ /^.{$len}ii/
         || $gap && $mutlinecp =~ /^.{$len}i\-i/ ) {
      --$mismatches;
      $adjustment += 9*$CtoTpenalty/10;
    }
    $len += 2;
    ++$len if $gap;
    $mutlinecp =~ s/^.{$len}//;
  }
  $adjustment = int($adjustment + 0.49); # rounding 100.5 down to 100 ; feels better;-)
  return ($adjustment,$mismatches);
}

sub calcambigadj {
  my $name = shift;
  my $seq = shift;
  $seq =~ s/\s//g;
  $seq =~ tr/a-z/A-Z/;
  my $match = $matscore{'A'}{'A'};
  $seq =~ tr/AT//d;
  my $gc = ($seq =~ tr/GC//d);
  if ($matscore{'G'}{'G'} != $match) {
    # earlier Xmatch.pl forbids C:C and G:G or A:A and T:T to have different scores                   
    $ambigadj{$name} += ($match - $matscore{'G'}{'G'})*$gc;
  }
  my $n= ($seq =~ tr/N//d);
  $ambigadj{$name} += ($match - $matscore{'N'}{'N'})*$n if $n;
  if ($seq) { # most seqs do not have other ambigs than N                                             
    my @ambigs = qw (R Y K M S W X B D H V Z);
    foreach my $a (@ambigs) {
      $ambigadj{$name} += ($match - $matscore{$a}{$a})*($seq =~ s/$a//g);
      # cannot use faster tr as it doesnt interpret variables                                         
      last unless $seq;
    }
  }
  if ($seq) {
    warn "\nNOTE: $name has irregular bases $seq\n\n";
  }
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
