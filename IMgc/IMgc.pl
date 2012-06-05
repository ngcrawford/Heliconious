#!/usr/bin/env perl
# 
# IMgc
#
# augustw@email.arizona.edu
# 
# this reads in a fasta file and makes an IM file out of it.
# the current strategy is to remove (a minimal # of) individuals that drive 4-gamete violations and
# infinite sites violations.

# bioperl requirement removed
# murray cox <m.p.cox@massey.ac.nz>
# 15 june 2011

use strict;
use Getopt::Long;

# these tell me if I should 1st try and weed out individuals, or if I should try splitting up the file
my ($WEED, $SPLIT, %segregatingSites, %fileGuts, @ids, @positions, @fastaIdsInOrder, 
  $INPUT_FILENAME, $SPLIT_CUTOFF, $help, $FASTA, $INDIVIDUAL_WEIGHT, $OMIT_RECORD_NUMBER);

$WEED = $SPLIT = 0;
# this is the weight factor for individuals and for sites 
# (ie, removing an individual is considered as bad as removing a segregating site it == 1)
# IE, this is the inclusiveness score parameter alpha
$INDIVIDUAL_WEIGHT = 1;

# I do the after-analysis on the top 10 fits
my $numberOfCandidatesToKeep = 10;
# the after analysis counts the literal number of sites that remain variable


# if you want to include a sequence, but not include them in detecting where the 4-gamete violations are, then
# flag this to the record number in the fasta(s)
$OMIT_RECORD_NUMBER = -1;
my $OMITTED_RECORD;
my $OMITTED_ID;
my $PRESERVED_RECORD_NUMBER = -1;
my $PRESERVED_ID;
# if your data has Ns in it, then this will optimize for including the 
# maximal number of sites where no individual contains an N
# (IM ignores the whole site if 1+ individuals have an N)
my $OPTIMIZE_FOR_NO_DATA = 0;

GetOptions("weed" => \$WEED,
		"fasta" => \$FASTA,
		"help" => \$help,
		"i=f" => \$INDIVIDUAL_WEIGHT,
		"o=i" => \$OMIT_RECORD_NUMBER,
		"p=i" => \$PRESERVED_RECORD_NUMBER,
		'n' => \$OPTIMIZE_FOR_NO_DATA,
		"split" => \$SPLIT);

my $usage = 
"Doh! Usage is: $0 name_of_fasta_file(s)
		(-w) [weeds out individuals; maintains original sequence length (fast)]
		(-s) [weeds out both individuals and sites; makes the largest block possible; (slow)]
		(-i 1.1 [weights keeping individuals. Default = 1. Weights keeping individuals vs keeping sites]
		(-o 0 [omits the nth id from your fasta file when optimizing the region, but tacks them back on in the output]
		(-p 0 [preserves the pth id from your fasta file no matter what]
		(-f) [prints output in fasta format (not IM)]\n";

if ($help) {
  die $usage;
}

if ($WEED && $SPLIT) {
  die "DoH! It's weed or split, not weed and split!\n" , $usage;
} elsif (!$WEED && ! $SPLIT) {
  $SPLIT = 1; # default
}

if ($PRESERVED_RECORD_NUMBER == $OMIT_RECORD_NUMBER &&
	$PRESERVED_RECORD_NUMBER != -1) {
  die "Doh! You can preserve record # $PRESERVED_RECORD_NUMBER or omit them; not both\n" , $usage;
} 

my $MAX_ID_SIZE = 9;

my $N_WARNING = 0;
my $N_WARNING_MESSAGE = 
"Just so you know. IM will read files with Ns, ?s and -s,
but it will not look at segregating sites with them\n\n";

if (!@ARGV || ! -e $ARGV[0]) {
  die $usage;
}


if (@ARGV == 1) {
  $INPUT_FILENAME = $ARGV[0];
} else {
# since this was originally designed for IM, 
# this allows the pairs of populations to be separate files
# (otherwise you have to make files out of all the pair-wise combinations of files
# which gets to be a bigggg and obnoxious numer).
# EG, you can run $0 with two fasta files (or 3 or 4...), as long as they align to each other
  $INPUT_FILENAME = &getCompositeFilename();

  if ($FASTA) {
    $INPUT_FILENAME .= '.fasta';
  }

}


# if the above is undef, then it's the number of fastaRecords/2

while (@ARGV) { 
  &loadFasta();
}


@ids = sort {$a cmp $b} keys %fileGuts; # the fasta record IDs from the input file

if (!&findSegSites()) {
  print STDERR "Sorry. No seg sites were found in your fasta file.\n";
  exit(0);
}

open OUTPUT, ">$INPUT_FILENAME.out";


if ($WEED) {
  &weedDatFile();
} elsif ($SPLIT) {  
 # &splitItTwo();
  &splitItThree();
} 

close OUTPUT;

if ($FASTA) {
  print "... whoo hoo! Fasta file made Aokay\n";
} else {
  print "... whoo hoo! IM file made Aokay\n";
}


sub splitItThree() {

  if (@positions < 2) {
      my %tempHash = ();
      print "There were < 2 segregating sites-- I'm just going to print the alignment block.\n";
      &printAlignmentBlock(0, undef,\%tempHash);      
      return;
  }
    
  print "There are " , scalar(@ids) , " ids and " , scalar(@positions) , " seg sites\n";
  # run the detection algorthm twice; once regularly, and once
  # treating the Ns as 4-gamete violations
  my $scoreRef = &splitItThreeHelper(0);
  if ($OPTIMIZE_FOR_NO_DATA) {
    my $scoreRefNsAsViolators = &splitItThreeHelper(1);
  # merge the arrays
    push(@$scoreRef, @$scoreRefNsAsViolators);
  }
  
  my $index = &findBestFit($scoreRef);
    # do a quick and dirty attempt at reinstating omitted individuals
  if ($index < @$scoreRef) {
    &relax($scoreRef->[$index]->[0], $scoreRef->[$index]->[1], $scoreRef->[$index]->[3], 0);
  } else {
    &relax($scoreRef->[$index]->[0], $scoreRef->[$index]->[1], $scoreRef->[$index]->[3], 1);
  }

  &printAlignmentBlock($scoreRef->[$index]->[0], $scoreRef->[$index]->[1], $scoreRef->[$index]->[3]);
}

sub relax($$$) {
  my ($startIndex, $endIndex, $perps, $nsAsViolators) = @_;
  my ($i, @keys, $testRemove);
  
  @keys = keys %$perps;
  
  IDS: while (@keys) {
    $testRemove = shift @keys; # grab a random subject
    delete $perps->{$testRemove}; # and try removing them
    
    $startIndex = $_[0];
    
    while ($startIndex < $endIndex)  { # look for 4-gamete violations, w/ this individual excluded
      for ($i = $startIndex + 1; $i < $endIndex; $i++) {
        if (ref &test4Gametes($startIndex, $i, $perps, $nsAsViolators)) { # four gamete violation
	  $perps->{$testRemove} = 1; # re-add the individual back
	  next IDS;
	}
      }
      $startIndex++;
    }
    print("Reinstating $testRemove!\n");
  }

}


sub splitItThreeHelper($) {
  my $nsAsViolators = shift;  
  my %fourGametePerpetrators = ();
  my ($startIndex, $endIndex, %candidates, $id, $key, @keys, $score, @scores, $i, $j);

  my $numberOfIndividuals = scalar @ids;

  @scores = %candidates = ();
 # test4Gametes
 # finds the four-gamete violations for all pairs of sites 
  for ($startIndex = 0; $startIndex < @positions; $startIndex++) {
    for ($endIndex = $startIndex + 1; $endIndex < @positions; $endIndex++) {
      $key = "$startIndex,$endIndex";
      $fourGametePerpetrators{$key} = &test4Gametes($startIndex, $endIndex, undef, $nsAsViolators);
      if ($fourGametePerpetrators{$key}) {
        foreach $id (@{$fourGametePerpetrators{$key}}) {
	  $candidates{$id} = 1;
	}
      }
    }
  }
    
  $score = $numberOfIndividuals - (scalar keys %candidates);
  
  # if there's no recombination, just return the whole block
  if ($score == $numberOfIndividuals) {
    $score **= $INDIVIDUAL_WEIGHT;
    $score *= @positions;
    push(@scores, [0, scalar @positions, $score, \%candidates]);
    return \@scores;
  }
 # recalculate the inclusiveness score 
  $score **= $INDIVIDUAL_WEIGHT;
  $score *= @positions;
  @scores = ();
  # the score for keeping the entire file, and just removing individuals
  push(@scores, [0, scalar @positions, $score, \%candidates]);

  my $bestScoreSoFar = $score;
  my $sitesToExclude = 1;
  my $maxPossibleScore;
  $startIndex = 0;
  # this progressively shrinks the alignment block
  # it finds the number of individuals needed to be removed for each 
  # continuous sub-alignment. It calculates a score for each subalignment
  # and after each shrinking, it looks at the next smallest size, and calculates
  # whether or not it's possible for the score to improve. If the score cannot 
  # get bigger, then the locally best scoring block is the globally best scoring block
  # and I break out of the loop
  
  INFINITE: while (1) {

    for ($startIndex = 0; $startIndex <=$sitesToExclude; $startIndex++) {
      $endIndex =  scalar @positions - $sitesToExclude + $startIndex;
      my %theseCandidates = ();
      # for a subsection of the alignment block, find the individuals you need to
      # remove in order to have a non-recombining fragment
      for ($i = $startIndex; $i < $endIndex; $i++) {
        for ($j = $i + 1; $j < $endIndex; $j++) {
	  $key = "$i,$j";
	  if ($fourGametePerpetrators{$key}) {
            foreach $id (@{$fourGametePerpetrators{$key}}) {
	      $theseCandidates{$id} = 1;
	    }  
	  }
	 
        }
      }

      @keys = keys %theseCandidates;
      $score = (($numberOfIndividuals - scalar(@keys))**$INDIVIDUAL_WEIGHT) * ($endIndex - $startIndex);
      push(@scores, [$startIndex, $endIndex, $score, \%theseCandidates]);
      if ($score > $bestScoreSoFar) {
        $bestScoreSoFar = $score;
      }
    }

    # assume that if you make the block 1-site smaller then you remove no individuals;
    # calculate the score for this possibility. If the current max is >= to this, then
    # exit the loop
    $maxPossibleScore = ($numberOfIndividuals ** $INDIVIDUAL_WEIGHT) * ($endIndex - $startIndex - 1);
    if ($maxPossibleScore < $bestScoreSoFar) {
      last INFINITE;
    }
    $sitesToExclude++;
  }
  
  return \@scores;
}


sub splitItTwo() {

  my $maximum = @ids * @positions;
  
  print "There are " , scalar(@ids) , " ids and " , scalar(@positions) , " seg sites\n";
  
  if (@positions < 2) {
      my %tempHash = ();
      print "There were < 2 segregating sites-- I'm just going to print the alignment block.\n";
      &printAlignmentBlock(0, undef,\%tempHash);      
      return;
  }
  
  my ($alignmentBlocks, $startIndex, $endIndex, $tolerance, 
  	$perpsHashRef, $comp1, $comp2, $sitesToExclude);
  
  $sitesToExclude = $tolerance = $startIndex = 0;
  
  my @candidates = ();

  print "I'm going to look for the " , ($numberOfCandidatesToKeep * 2) , " best fits possible, and then pick the prettiest of them all\n";

# this walks the parameter space, attempting to maximize the amount of data in the
# recombining blocks found
  INFINITE: while (1) {
  # perps are those that must be removed to have an alignment block
   for ($startIndex = 0; $startIndex <=$sitesToExclude; $startIndex++) {
     $endIndex =  scalar @positions - $sitesToExclude + $startIndex;
#      print "$startIndex $endIndex\n";
     ($alignmentBlocks, $perpsHashRef) = &fourGameteThis($startIndex, $endIndex, $tolerance, 0);
     
     if ($alignmentBlocks == 1) {
       push(@candidates, [$startIndex, $endIndex, $tolerance, $perpsHashRef]);
       if (@candidates >= $numberOfCandidatesToKeep) {
         last INFINITE;
       }
     }
   }
# compare the cost of removing and indivudal from the alignment to the cost of
# removing an additional site.
    $comp1 = (((scalar(@ids) / $INDIVIDUAL_WEIGHT) - $tolerance - 1) * @positions) / $maximum;
    $comp2 = (scalar(@ids) * (@positions - $sitesToExclude - 1)) / $maximum;
    # if tolerating removing more individuals will preserve more of the file, then try it
    if ($comp1 >= $comp2) {
      $tolerance++;
    } else { # otherwise, move another site
      $sitesToExclude++;
    }
  }
  
  $sitesToExclude = $tolerance = 0;
  $numberOfCandidatesToKeep *= 2;
  
  # IM treats sites w/ Ns as non-existent. This attempts to account for this fact
  # let it rewalk the parameter space, but treat Ns as 4-gamete violations, and see what you find
  INFINITE: while (1) {
  # perps are those that must be removed to have an alignment block
   for ($startIndex = 0; $startIndex <=$sitesToExclude; $startIndex++) {
     $endIndex =  scalar @positions - $sitesToExclude + $startIndex;
#      print "$startIndex $endIndex\n";
     ($alignmentBlocks, $perpsHashRef) = &fourGameteThis($startIndex, $endIndex, $tolerance, 0, 1);
     
     if ($alignmentBlocks == 1) {
       push(@candidates, [$startIndex, $endIndex, $tolerance, $perpsHashRef]);
       if (@candidates >= $numberOfCandidatesToKeep) {
         last INFINITE;
       }
     }
   }
    
    $comp1 = (((scalar(@ids) / $INDIVIDUAL_WEIGHT) - $tolerance - 1) * @positions) / $maximum;
    $comp2 = (scalar(@ids) * (@positions - $sitesToExclude - 1)) / $maximum;
    # if tolerating removing more individuals will preserve more of the file, then try it
    if ($comp1 >= $comp2) {
      $tolerance++;
    } else { # otherwise, move another site
      $sitesToExclude++;
    }
  }
  
  my $bestFitIndex = &findBestFit(\@candidates);
  print "The prettiest is from position $candidates[$bestFitIndex][0] to $candidates[$bestFitIndex][1] ignoring up to $candidates[$bestFitIndex][2] individuals\n";
  
  &printAlignmentBlock($candidates[$bestFitIndex][0], $candidates[$bestFitIndex][1], $candidates[$bestFitIndex][3]);
}

# this takes in the array of candidates from splitItTwo(), and finds the best one given the actual data preserved
sub findBestFit($) {
  my $arrayRef = shift;
  my ($numberOfUsedSites, $numberOfIds, $bestIndex, $productMax, $index, $idSize);
  $bestIndex = $numberOfUsedSites = $numberOfIds = $productMax = 0;
  my $i = 0;
  foreach my $candidate (@$arrayRef) {
    for ($index = $candidate->[0]; $index < $candidate->[1]; $index++) {
      $numberOfUsedSites += &isVariable($index, $candidate->[3]);
    }
    
    $idSize = scalar(@ids) - scalar(keys(%{$candidate->[3]}));
    print "$i: $numberOfUsedSites sites for $idSize individuals\n";

    # recycle the variable name and 
    # calculate the inclusiveness score
    $idSize **= $INDIVIDUAL_WEIGHT;
    $idSize *= @positions;
    if ($idSize > $productMax) {
      $productMax = $idSize;
      $bestIndex = $i;
    }
    $numberOfUsedSites = 0;
    $i++;
  }
  print "Index $bestIndex picked\n";
  return $bestIndex;
}

# given an index in the seg-sites, and a hash of offenders, tell me if IM can use the site 
# IM can't use indels, and can't use sites w/ Ns, and won't use sites who are only variable in the offenders
sub isVariable($$) {
  my ($index, $offenders) = @_;
  my $base1 = "";
  my $isVariable = 0;
  
  foreach my $id (@ids) {
    if (exists $offenders->{$id}) {
      next;
    }
    
    if ($segregatingSites{$id}[$index] =~ m/[N-]/ ) {
      # if there's no data, then effectively the site is not variable
      if ($OPTIMIZE_FOR_NO_DATA) { 
        return 0;
      }
    } elsif ($base1 eq "") {
      $base1 = $segregatingSites{$id}[$index];
    } else {
      $isVariable = 1;
    }
  }
  
  return $isVariable;
}


# this is what minimizes the recombination events...
# for the weed pragma, it necessarily keeps the entire fasta file, and converts the whole thing 
# to n IM-compatible alignment blocks
# The number of alignment blocks is inversely related to the number of individuals in the fasta file-- thus
# this gives the users a sliding set of stringencies; eg, it will query the user if they want, say, 5 alignment blocks
# if 3 individuals are removed, 4 alignment blocks if 8 individuals are removed... etc.
#
# note that the target # of individuals can be different between different alignment blocks
sub weedDatFile() {
  if (@positions < 2) {
      my %tempHash = ();
      print "There were < 2 segregating sites-- I'm just going to print the alignment block.\n";
      &printAlignmentBlock(0, undef,\%tempHash);      
      return;
  }
  
  my (@blocks, $returnValue);
  @blocks = ();
  print "Preparing to look for singleton perpetrators to remove\n";

  for (my $i = 0; $i <= $SPLIT_CUTOFF; $i++) {
    $returnValue = &fourGameteThis(0, scalar(@positions), $i, 0);
    
    if (!@blocks ||
    	$returnValue < $blocks[-1][1]) {
      push(@blocks, [$i, $returnValue]);
    }
    
    # if I can get it down to one block, then that's as good as it can get.
    if ($returnValue < 2) {
      last;
    }
     
  }
  # get the user-feedback, for what they consider to be the optimal tradeoff between removing individuals, and splitting the alignment blocks
  my $userFeedback = &whichIndex(\@blocks);

# use the user-defined tolerance, and tell it to print  
  &fourGameteThis(0, scalar(@positions), $blocks[$userFeedback][0] , 1);
}


sub fourGameteThis($$$$$) {
 # the start and stop index to look at for 4-gamete violations,
 # the tolerance for the number of individuals that can be removed,
 # and whether or not to print the results to file
  my ($segSiteStart, $segSiteEnd, $tolerance, $print, $nsAsViolators) = @_;
  # who has caused the 4-gamete violations, the 2 indeces to be compared for the 4-gamete violation test
  # the scalar of the %
  # and the return value from the test4gamete sub
  my ($index1, $index2, $violations, $returnValue, $id);
  # @returnValue = (isFourGameteViolation (boolean)
  #                  isindex1lost
  #                  isIndex2lost
  # OR
  # @returnValue = (arrayOfAdditional4GameteViolators
  #                  isindex1lost
  #                  isIndex2lost
  
  $index1 = $segSiteStart;
  $index2 = $index1 + 1;
  
  if (!defined $nsAsViolators) {
    $nsAsViolators = 0;
  }
  
  my %violators = ();
  
  my $alignmentBlocks = 1;
  
  # while loop used to do the pair-wise comparisons
  while ($index2 < $segSiteEnd) {
    # note-- %violators is not modified by test4gametes
 #   print "$index1 vs $index2\n";   

    $returnValue = &test4Gametes($index1, $index2, \%violators, $nsAsViolators);
    
    # if it's a reference, then it's a reference to the individual perps
    if (ref $returnValue) {
      $violations += scalar(@$returnValue);
      
      if ($violations > $tolerance) {
        
        if ($SPLIT) {
        # if I'm trying to get the biggest block, then this set of tested indeces didn't work
          return (2, \%violators);
        } elsif ($print && $WEED) { #if I'm trying to maximize the inclusiveness of the input file, print the current block
          &printAlignmentBlock($segSiteStart, $index2, \%violators);
        }
        
        $alignmentBlocks++;
        $violations = 0;
        %violators = ();
        # start trying to grab another alignment block
        $segSiteStart = $index2;
        $index1 = $index2;
        $index2++;
        next;
      } else {
      # if we are still below the quota, keep trying to get more sites, and record who the perps are
        foreach my $id (@$returnValue) {
          $violators{$id} = "";
        }        
      }
    }

    $index1++;

    if ($index1 >= $index2) {
      $index2++;
      $index1 = $segSiteStart;      
    }
    
  } # end while loop
  
  if ($print) {
    &printAlignmentBlock($segSiteStart, $segSiteEnd, \%violators);
  } 
  
  if ($SPLIT) {
    return ($alignmentBlocks, \%violators);
  }
  
  return $alignmentBlocks;
}

# this does the 4-gamete test on 2 sites for all individuals.
# 
sub test4Gametes($$$$) {
  my $index1 = shift;
  my $index2 = shift;
  my $violatorsRef = shift;
# lets' me optionally consider individuals w/ Ns as four-gamete violations
  my $noDataFlag = shift;
  my @Ns;
  
  if (defined $noDataFlag && $noDataFlag) {
    @Ns = (); # initialize the array of No data
  } else {
    $noDataFlag = 0;
  }

  
  my $id;
  my %characterStates = ();

  foreach $id (@ids) {
    
    # if the individual is already trimmed from the outputfile to be, then ignore 'em
    if (defined $violatorsRef && exists $violatorsRef->{$id}) {
      next;
    }
    
    # ignore no data
    if ($segregatingSites{$id}[$index1] =~ m/[N?]/ || 
          $segregatingSites{$id}[$index2] =~ m/[N?]/) { 
      if ($noDataFlag && ! grep {$_ eq $id} @Ns) { # grab the unique IDs
        push(@Ns, $id);
      }	  
      next;
    }
            
    # associate the individuals to the characterstates of the two snps being compared
    if (!exists $characterStates{$segregatingSites{$id}[$index1] . "," . $segregatingSites{$id}[$index2]}) {
      $characterStates{$segregatingSites{$id}[$index1] . "," . $segregatingSites{$id}[$index2]} = [];
    }
    
    push(@{$characterStates{$segregatingSites{$id}[$index1] . "," . $segregatingSites{$id}[$index2]}}
              , $id); 
  }
  # sort by frequency
  my @keys = sort { scalar(@{$characterStates{$a}}) <=> 
  				scalar(@{$characterStates{$b}})
					} keys %characterStates;
   
  
  if (@keys < 4) {
  # return the individuals w/ No Data
    if ($noDataFlag && @Ns) {
      if (grep {$_ eq $PRESERVED_ID} @Ns) {
        return 0;
      }
      return \@Ns;
    }
    return 0;
  } elsif (@keys > 4) {
    die "Dooh! Fatal 4-gamete error! @keys\n";
  }
#  print " @{$characterStates{$keys[0]}}\n";

  
  if ($PRESERVED_RECORD_NUMBER == -1) {
    if ($noDataFlag && @Ns) {
      push(@{$characterStates{$keys[0]}}, @Ns);
    }
    return \@{$characterStates{$keys[0]}};
  }
  
  foreach my $key (@keys) {
  # return the least-common grouping that does not contain the preserved individual
    if (!grep {$_ eq $PRESERVED_ID} @{$characterStates{$key}}) {
      if ($noDataFlag &&
      		! grep {$_ eq $PRESERVED_ID} @Ns) {
     	push(@{$characterStates{$key}}, @Ns); 
      }
      return \@{$characterStates{$key}};
    }
  }
  die "Baad!\n";
}

# finds the segregating sites from %fileGuts 
sub findSegSites() {
  @positions = (); # the substring positions of the segregating sites.
  my $id;
  my $length = length($fileGuts{$ids[0]});
  my $substringPosition = 0;
  my $isSegregating;
  my $previousBase;
  my $base;
  my %characterStates;
  my @characterKeys;
  my $sizeOfEvent; # how big the  polymorphism is (1 = SNP, indels can be multibase)
  
  while ($substringPosition < $length) {
    $previousBase = $base = "";
    $isSegregating = 0;
    $sizeOfEvent = 1; # assume SNP (initially)
    
    foreach $id (@ids) {
       $base = substr($fileGuts{$id}, $substringPosition, 1);
       
       # tests for indels
       if ($base eq '-') {
         my $tempBase = $base;
	 my $thisEventSize = 0;
	 # greedily grab as many bases of the indel as possible
	 while (++$thisEventSize + $substringPosition < $length && $tempBase eq '-') {
	   $tempBase = substr($fileGuts{$id}, $substringPosition + $thisEventSize, 1);
	 } 
	 $thisEventSize--;
	 if ($thisEventSize > $sizeOfEvent) {
	   $sizeOfEvent = $thisEventSize;
	 }
       }
       
       if ($base eq $previousBase) {
         next;
       } elsif ($base eq 'N') { # ignore no data
         next;
       } elsif ($previousBase eq "") {
         $previousBase = $base;
       } else {
         $isSegregating = 1;
	 last;
       }
               
    }
    
    if (!$isSegregating) {
      $substringPosition+=$sizeOfEvent;
      next;
    }
# if it's a segregating site, record it

    print $substringPosition , " is a segregating site position\n";
    push(@positions, $substringPosition);
    
    %characterStates = ();

    foreach $id (@ids) {
       $base = substr($fileGuts{$id}, $substringPosition, $sizeOfEvent);
# change any site that contains an N to an N (for simplicity's sake)
       if ($base =~ m/[N?]/) {
 
	 if (!$N_WARNING) {
           print "$base found \n" , $N_WARNING_MESSAGE;
      	   $N_WARNING = 1;
         }
	 # if it's a mixture of Ns and GATC-, then change the whole thing to Ns
	 if ($base !~ m/^[N?]+$/) {
	   warn $id , " has No data interspersed with real data within an indel at site $substringPosition\n" ,
	   " I'm just going to call the whole event an N for this site.\n";
	   $fileGuts{$id} = substr($fileGuts{$id}, 0, $substringPosition) . ( 'N' x $sizeOfEvent ) .
	  	substr($fileGuts{$id}, $substringPosition + $sizeOfEvent);
	 }
         $base = 'N'; 
       } elsif (!exists $characterStates{$base}) {
         $characterStates{$base} = 1; # keep track of the quantity of bases
       } else {
         $characterStates{$base}++;
       }

       # segSites is a hash of arrays
       if (!exists $segregatingSites{$id}) {
         $segregatingSites{$id} = [$base];
       } else {
         push(@{$segregatingSites{$id}}, $base);
       }
    }
    # sort by basecall frequency
    @characterKeys = sort { $characterStates{$b} <=> $characterStates{$a} } keys %characterStates;
    
        
    if (@characterKeys < 3) {
     $substringPosition+=$sizeOfEvent;
      next;
    }
    print "Infinite sites violation found. I'll try and fix it.\n";
    # fix infinite sites violations...

    for (my $i = 2; $i < @characterKeys; $i++) {

      foreach $id (@ids) {
        if ($segregatingSites{$id}[-1] eq $characterKeys[$i]) {
	# if there's 3+ character states, change the lowest frequency one(s) to Ns in the seg-sites table
	# and in the body of the sequence data
	print "Removing one base call for $id at site $substringPosition because of an infinite sites violation\n";
	  $segregatingSites{$id}[-1] = 'N';
	  $fileGuts{$id} = substr($fileGuts{$id}, 0, $substringPosition) . ( 'N' x $sizeOfEvent ) .
	  	substr($fileGuts{$id}, $substringPosition + $sizeOfEvent);
	}
      }
    }
   $substringPosition+=$sizeOfEvent;
  }
  return scalar(@positions);
}


# loads the contents of the fasta file into a hash; id = key, sequence = value
sub loadFasta() {
  
  print "Loading $ARGV[0]\n";
  my $input_file = shift(@ARGV);
  
  if( -e $input_file ){
  	open( FASTA, "<$input_file" ) or die "error: file $input_file cannot be opened: $!\n";
  }else{
    die "error: file $input_file does not exist\n";
  }
  
  my $id = "";
  my $sequence = "";
  my $length = "";
  my $counter = 0;
  $OMITTED_ID = $OMITTED_RECORD = "";
  $PRESERVED_ID = "";
  
  # check format
  my $file_format = &get_format(*FASTA);
  if( $file_format ne "fasta" ){
    die "error: file $input_file is not in FASTA format\n";
  }
  
  until( eof(FASTA) ){
  
	my @sequence_returns = &get_next(*FASTA, $file_format);
	
	$id = $sequence_returns[0];
	$sequence = $sequence_returns[1];
	$length = length($sequence);
	
    if ($counter == $OMIT_RECORD_NUMBER) {
      $OMITTED_RECORD = uc($sequence);
      $OMITTED_ID = $id;
      print "Omitting $id\n";
      $counter++;
      next;
    } elsif ($counter == $PRESERVED_RECORD_NUMBER) {
      print "Preserving $id\n";
      $PRESERVED_ID = $id;
    }
    
    if (exists $fileGuts{$id}) {
      $id .= '_B';
    }
        
    if (!exists $fileGuts{$id}) {
      $fileGuts{$id} = uc($sequence);
      push(@fastaIdsInOrder, $id );    # lets me keep the original order of the IDs in the file 
      
      if (length $fileGuts{$fastaIdsInOrder[0]} != $length ) {
        print STDERR "Doh! Your fasta files must be aligned! The error occurred in " , $id , "\n";
        close FASTA or die "error: file $input_file cannot be closed: $!\n";
        die $usage;
      }
      
    } else {
      print STDERR "Doh! I need unique ids! ", $id , " is repeated\n";
      close FASTA or die "error: file $input_file cannot be closed: $!\n";
      die $usage;
    }
    $counter++;
  }

  close FASTA or die "error: file $input_file cannot be closed: $!\n";
  
  # if the user selected an index that is out of bounds, then...
  if ($OMIT_RECORD_NUMBER >= 0 && ! $OMITTED_ID) {
    die "Doh! You tried to omit record number $OMIT_RECORD_NUMBER, but that is not a valid index in this file\n";
  } elsif (!$counter) {
    die "Doh! Your fasta file has no sequences in it\n";
  } elsif ($counter == 1 && $OMITTED_ID) {
    die "Your fasta file has just one sequence in it, and it's the omitted record. That's just mean\n";
  }
  
  $SPLIT_CUTOFF = int(@fastaIdsInOrder/2);

  return 1;
}

sub printAlignmentBlock($$$) {
# $startFrom, $index2, \%offenders
  my ($startFrom, $endAt, $offendersRef) = @_;
  my $id;
  my $numberOfSpaces;
  # if start from is from 0, we want the substring to be from 0 (the sequence preceding the 1st segregating site too)
  # otherwise, we want the substring to start at the nth segregating site, and go thru to the base before the nth + x seg site

  if (defined $endAt && $endAt < @positions) {
    print "Seg sites: $startFrom to $endAt $positions[$endAt]\n";
  } else {
    print "Seg sites: $startFrom to end\n";
  }
  
  if ($startFrom) {
    $startFrom = $positions[$startFrom];
  }
  
    # re-include the previously excluded individual...
  if ($OMITTED_RECORD) {
# print the record identifier
    if ($FASTA) {
      print OUTPUT '>' , $OMITTED_ID , "\n";    
    } else {
        if (length($OMITTED_ID) > $MAX_ID_SIZE) {   
          print OUTPUT substr($OMITTED_ID, 0, $MAX_ID_SIZE) , " ";
        } else {
          $numberOfSpaces = $MAX_ID_SIZE - length($OMITTED_ID) + 1;
	  print OUTPUT $OMITTED_ID , (" " x $numberOfSpaces);
        }
    }
    # print the data
    if (defined $endAt && $endAt < @positions) {
      print OUTPUT substr($OMITTED_RECORD, $startFrom, $positions[$endAt] - $startFrom) , "\n";
    } else {
      print OUTPUT substr($OMITTED_RECORD, $startFrom) , "\n";
    }
  }
  
  
  foreach $id (@fastaIdsInOrder) {
  # format the ID printing portion
    if (!exists $offendersRef->{$id}) {
      # prints the header fasta-style
      if ($FASTA) {
        print OUTPUT '>' , $id , "\n";
      } else {
    
        if (length($id) > $MAX_ID_SIZE) {   
          print OUTPUT substr($id, 0, $MAX_ID_SIZE) , " ";
        } else {
          $numberOfSpaces = $MAX_ID_SIZE - length($id) + 1;
	  print OUTPUT $id , (" " x $numberOfSpaces);
        }
      }
# if the ending index is defined, and it's less than the # of seg-sites we have, then use it as a substring position      
      if (defined $endAt && $endAt < @positions) {
        print OUTPUT substr($fileGuts{$id}, $startFrom, $positions[$endAt] - $startFrom) , "\n";
      } else {
        print OUTPUT substr($fileGuts{$id}, $startFrom) , "\n";
      }
      
    } else {
      print "Excluding $id\n";
    }
  
  }
  

  
 
  if (defined $endAt && $endAt < @positions) {    
    print "AKA, substring position $startFrom to $positions[$endAt] \n"; 
  } else {
    print "AKA, substring position : $startFrom to end\n";
  }
  
  print OUTPUT "\n";

}

# parses through the argv and tries to make a good filename to describe all files therein.
sub getCompositeFilename() {
  my @argv_clone = @ARGV;
  
  my $arg;
  my @components;
  my @split;
  my $value;
  # makes an array of the non-repetitive bits of information between the N files, and makes an in-order string out of it
  foreach $arg (@argv_clone) {
    $arg =~ s/\.[^.]+$//;  # ditch the file-extension
    @split = split /[-_]+/, $arg;
    
    foreach $value (@split) {
      if (!grep {uc($value) eq uc($_)} @components) {
        push(@components, $value);
      }
    }
        
  }
  
  return join '_', @components;
  
}


sub whichIndex($) {
  my $blockRef = shift;
  my $feedback;

  print "The following minima are available\nPlease select by index\n";
  
  
  my $i = 0;
  
  foreach my $arrayRef (@$blockRef) {
    $i++;
   

      print $i , ": If you remove $arrayRef->[0] individual";
    
      if ($arrayRef->[0] != 1) {
        print 's';
      }
      print " you will get $arrayRef->[1] alignment block";
   
      if ($arrayRef->[1] != 1) {
        print 's';
      }
    
      print "\n";
    
  }
 
  print ":";
  
  while (!defined $feedback) {
    $feedback = <STDIN>;
    chomp($feedback);
    # if the user-input is sensical
    if ($feedback =~ m/^\d+$/ && $feedback != 0 && $feedback <= @$blockRef) {
      $feedback--;
      return $feedback;
    } else {
      print "Eh? I don't get it.\n";
      $feedback = undef;
    }
    
  }
  
}

# determine FASTA/FASTQ format
# usage: &get_format(*FILEHANDLE)
# return: "fasta" or "fastq"
#         0 on failure
sub get_format(*){
	
	# set function variables
	local *FILEHANDLE = shift;
	
	# retrieve file position
	my $position = tell FILEHANDLE;
	
	# retrieve first line
	seek(FILEHANDLE, 0, 0);
	my $first_line = <FILEHANDLE>;
	
	# retrieve first character
	my $first_character = substr($first_line, 0, 1);
	
	# reset filehandle
	seek(FILEHANDLE, $position, 0);
	
	# return format
	if( $first_character eq ">" ){
		return "fasta";
	}elsif( $first_character eq "@" ){
		return "fastq";
	}else{
		return 0;
	}
}

# retrieve next FASTA/FASTQ entry
# usage: &get_next(*FILEHANDLE, $format)
# return: @array
#         $array[0] = header, $array[1] = sequence, $array[2] = quality
#         0 on failure
sub get_next(*$){
	
	# set function variables
	local *FILEHANDLE = shift;
	my $format = shift;
	
	my @data = ("", "", "");
	
	if( $format eq "fasta" ){
		
		# retrieve first line
		my $first_line = <FILEHANDLE>;
		chomp($first_line);
		
		# retrieve file position
		my $position = tell FILEHANDLE;
		
		# retrieve header
		$data[0] = substr($first_line, 1, length($first_line)-1);
		
		# retrieve sequence
		until( eof(FILEHANDLE) ){
			
			# retrieve line
			my $line = <FILEHANDLE>;
			chomp($line);
			
			# retrieve first character
			my $first_character = substr($line, 0, 1);
			
			# step through multiline fasta
			if( $first_character eq ">" ){
				seek(FILEHANDLE, $position, 0);
				last;
			}else{
				$data[1] .= $line;
				$position = tell FILEHANDLE;
			}
		}
	}elsif( $format eq "fastq" ){
		
		# retrieve first line
		my $first_line = <FILEHANDLE>;
		chomp($first_line);
		
		# retrieve header
		$data[0] = substr($first_line, 1, length($first_line)-1);
		
		# retrieve sequence
		my $second_line = <FILEHANDLE>;
		chomp($second_line);
		$data[1] = $second_line;

		# ignore next header line
		my $third_line = <FILEHANDLE>;

		# retrieve quality
		my $fourth_line = <FILEHANDLE>;
		chomp($fourth_line);
		$data[2] = $fourth_line;
	}else{
		return 0;
	}
	
	# return data
	return @data;
}
