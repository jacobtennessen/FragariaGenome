#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_t $opt_s);

# Usage
my $usage = "
AssessFragariaTrees.pl
Reads a Newick format tree and records position of target taxon

Copyright (C) 2019 by Jacob A Tennessen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: perl AssessFragariaTrees.pl options
 required:
  -t	(path to) tree file directory
  -s  coma-delimited list of target taxon names
";

#############

getopts('t:s:');
die $usage unless ($opt_t);
die $usage unless ($opt_s);

my ($directory, $samplelist);

$directory = $opt_t if $opt_t;

my @dirdata = split /\//, $directory;

$directory = join "/", @dirdata;

$samplelist = $opt_s if $opt_s;

my @samples = split ",", $samplelist;

my %colors = (
  "vesca" => "#e41a1c",
  "viridis" => "#4daf4a",
  "nipponica" => "#984ea3",
  "iinumae" => "#377eb8",
  "Camarosa_vesca" => "#f8c1c2",
  "Camarosa_viridis" => "#cce9cb",
  "Camarosa_nipponica" => "#e0c4e4",
  "Camarosa_iinumae" => "#bfd6ec",
  "Camarosa_iinumae,iinumae" => "#214e70",
  "Camarosa_nipponica,nipponica" => "#5e3067",
  "Camarosa_viridis,viridis" => "#2f6d2e",
  "Camarosa_vesca,vesca" => "#901010",
  "Camarosa_nipponica,Camarosa_viridis" => "#a65628",
  "Camarosa_iinumae,Camarosa_viridis,iinumae" => "#82b4c2",
  "Camarosa_iinumae,Camarosa_nipponica,iinumae" => "#8ca1ce",
  "Camarosa_iinumae,Camarosa_viridis" => "#c6e0dc",
  "Camarosa_iinumae,Camarosa_nipponica" => "#d0cde8",
);

my %chromlengths = (
  "Fvb1" => 24253023,
  "Fvb2"=> 29351230,
  "Fvb3" => 38324302,
  "Fvb4" => 33907851,
  "Fvb5" => 29430145,
  "Fvb6" => 39795230,
  "Fvb7" => 24229589,
);

my %targetdiploidcolors = (
  "vesca" => "#e41a1c",
  "bracteata" => "#e41a1c",
  "mandshurica" => "#f8c1c2",
  "bucharica" => "#f8c1c2",
  "nilgerrensis" => "#ff7f00",
  "viridis" => "#4daf4a",
  "nipponica" => "#984ea3",
  "iinumae" => "#377eb8",
);

my %targetdiploidcolorsfaded = (
  "vesca" => "#f8c1c2",
  "bracteata" => "#f8c1c2",
  "mandshurica" => "#f8c1c2",
  "bucharica" => "#f8c1c2",
  "nilgerrensis" => "#ffdab4",
  "viridis" => "#cce9cb",
  "nipponica" => "#e0c4e4",
  "iinumae" => "#bfd6ec",
);


my @origcolortaxa = sort (keys %colors);

my @origcolors;

my @targetdiploids = ("vesca","nilgerrensis","viridis","nipponica","iinumae");

foreach my $td (@targetdiploids) {
  push @origcolors, $targetdiploidcolors{$td};
  push @origcolors, $targetdiploidcolorsfaded{$td};
}

push @origcolors, "black";

my $spacer = 5000000;

my $cumul = 0;

my %position;

foreach my $c (sort keys %chromlengths) {
  $position{$c} = $cumul;
  $cumul += $chromlengths{$c};
  my $cmid = $position{$c} + ($cumul-$position{$c})/2;
  print "Chrom $c from $position{$c} to $cumul, mid is $cmid\n";
  $cumul += $spacer;
}

print "Spaced genome is $cumul in length\n";

my %fullsisters;

my $samplecount = 0;
  
my @percentout;

my $percentoutfile = "$directory/Percent_Sisters.txt";

my @percenttitle = ("Chrom","Subgenome");

foreach my $s (@origcolors) {
    push @percenttitle, $s;
}

my $percenttitle = join "\t", @percenttitle;

push @percentout, $percenttitle;

foreach my $sample (@samples) {
  
  print "\n#Examining $sample\n";
  
  my %localsisters;
  
  my %localdip;

  my @out;
  
  push @out, "Chrom\tStart\tEnd\tPosition\tSister\tColor";
  
  my @outdip;
  
  push @outdip, "Chrom\tStart\tEnd\tPosition\tClosestDiploid\tColor";
  
  my $outfile = "$directory/$sample"."_tree_windows.txt";
  
  my $dipoutfile = "$directory/$sample"."_tree_windows_closest_diploid.txt";
  
  my $windowcount = 0;
  
  foreach my $chrom (sort keys %chromlengths) {
    
    my %chromlocaldip;
    
    my $chromwindowcount = 0;
    
    my $treefile = "$directory/$chrom.masked-hets.raxml";
  
    if (open(IN, $treefile)) {
      while (<IN>) {
        my $line = $_;
        $line =~ s/\r|\n//g;
        my @metadata = split "\t", $line;
        if ($metadata[4] =~ /\(/) {
          my $goodsister;
          my $closestdip;
          my $dipcolor = "black";
          if ($metadata[4] =~ /$sample/) {
            my @data = split "", $line;
            my $incount = 0;
            my $pos = 0;
            my %HoAworkingpos;
            my %HoAsistersbyposition;
            my %HoAnestbyposition;
            my %clades;
            my %HoHsistersbysister;
            my %levelcount;
            foreach my $d (@data) {
              if ($d =~ /\(/) {
                $incount +=1;
                if (defined $levelcount{$incount}) {
                  $levelcount{$incount} +=1;
                } else {
                  $levelcount{$incount} = 1;
                }
                @{$HoAworkingpos{$incount}} = ();
                push @{$HoAworkingpos{$incount}}, $pos + 1;
              } elsif (($d =~ /,/)||($d =~ /\)/)) {
                my @cladedata = @data[$HoAworkingpos{$incount}[-1]..$pos];
                until ((!defined $cladedata[0])||(($cladedata[-1] !~ /\)/)&&($cladedata[-1] !~ /,/))) {
                  my $junk = pop @cladedata;
                }
                my $clade = join "", @cladedata;
                push @{$HoAsistersbyposition{$incount}}, $clade;
                for (my $i = 1; $i <= $incount; $i++) {
                  my $level = "$i\t$levelcount{$i}";
                  push @{$HoAnestbyposition{$level}}, $clade;
                }
                $clades{$clade} = $incount;
                if ($d =~ /,/) {
                  push @{$HoAworkingpos{$incount}}, $pos + 1;
                } elsif ($d =~ /\)/) {
                  my @seen;
                  foreach my $c (@{$HoAsistersbyposition{$incount}}) {
                    foreach my $s (@seen) {
                      $HoHsistersbysister{$c}{$s} = 1;
                      $HoHsistersbysister{$s}{$c} = 1;
                    }
                    push @seen, $c;
                  }
                  @{$HoAsistersbyposition{$incount}} = ();
                  $incount -=1;
                }
              }
              $pos +=1;
            }
            my $highestclade = 0;
            foreach my $level (keys %HoAnestbyposition) {
              my $seendiploid;
              my $seensample;
              foreach my $clade (@{$HoAnestbyposition{$level}}) {
                if ($clade =~ /$sample/) {
                  $seensample = 1;
                }
                if (defined $targetdiploidcolors{$clade}) {
                  if (defined $seendiploid) {
                    $seendiploid = "$seendiploid,$clade";
                  } else {
                    $seendiploid = $clade;
                  }
                }
                if ((defined $seendiploid)&&(defined $seensample)) {
                  my @leveldata = split "\t", $level;
                  if ($leveldata[0] > $highestclade) {
                    $highestclade = $leveldata[0];
                    $closestdip= $seendiploid;
                  }
                }
              }
            }
            foreach my $c (keys %clades) {
              if ($c =~ /^$sample:/) {
                my @sampledata = split ":", $c;
                my $shortest;
                foreach my $sister (keys %{$HoHsistersbysister{$c}}) {
                  my @fullsisterdata = split ":", $sister;
                  my @sisterdata;
                  foreach my $fsd (@fullsisterdata) {
                    my @fsd1 = split /\(/, $fsd;
                    if ((defined $fsd1[-1])&&(length($fsd1[-1]) >= 2)) {
                      my @fsd2 = split ",", $fsd1[-1];
                      if ((defined $fsd2[-1])&&(length($fsd2[-1]) >= 2)) {
                        my @fsd3 = split /\)/, $fsd2[-1];
                        if ((defined $fsd3[0])&&(length($fsd3[0]) >= 2)&&($fsd3[0] =~ /\D\D/)) {
                          push @sisterdata, $fsd3[0];
                        }
                      }
                    }
                  }
                  my $sisterline = join ",", sort (@sisterdata);
                  unless ((defined $shortest)&&($shortest < (scalar(@sisterdata)))) {
                    $shortest = (scalar(@sisterdata));
                    $goodsister = $sisterline;
                  }
                }
              }
            }
            unless (defined $goodsister) {
              print "Odd topology at $line\n";
              exit;
            }
            unless (defined $colors{$goodsister}) {
              $colors{$goodsister} = "black";
            }
            if ((defined $closestdip)&&(defined $targetdiploidcolors{$closestdip})) {
              if ($closestdip =~ /$goodsister/) {
                $dipcolor = $targetdiploidcolors{$closestdip};
              } else {
                $dipcolor = $targetdiploidcolorsfaded{$closestdip};
              }
            }
          } else {
            $goodsister = "NA";
            $colors{$goodsister} = "grey70";
            $closestdip = "NA";
            $dipcolor = "grey70";
          }
          my $end = $metadata[1] + $metadata[2];
          my $pos = $position{$metadata[0]} + $metadata[1] + $metadata[2]/2;
          push @out, "$metadata[0]\t$metadata[1]\t$end\t$pos\t$goodsister\t\"$colors{$goodsister}\"";
          push @outdip, "$metadata[0]\t$metadata[1]\t$end\t$pos\t$closestdip\t\"$dipcolor\"";
          if (defined $fullsisters{$goodsister}) {
            $fullsisters{$goodsister} +=1;
          } else {
            $fullsisters{$goodsister} = 1;
          }
          if (defined $localsisters{$goodsister}) {
            $localsisters{$goodsister} +=1;
          } else {
            $localsisters{$goodsister} = 1;
          }
          if (defined $localdip{$dipcolor}) {
            $localdip{$dipcolor} +=1;
          } else {
            $localdip{$dipcolor} = 1;
          }
          if (defined $chromlocaldip{$dipcolor}) {
            $chromlocaldip{$dipcolor} +=1;
          } else {
            $chromlocaldip{$dipcolor} = 1;
          }
          unless ($goodsister =~ /^NA$/) {
            $windowcount +=1;
            $chromwindowcount +=1;
          }
        }
      }
      close (IN);
    }
    
    my @percentline = ($chrom,$sample);
    
    foreach my $s (@origcolors) {
      if (defined $chromlocaldip{$s}) {
        my $cpercent = sprintf "%.4f", (100*$chromlocaldip{$s}/$chromwindowcount);
        push @percentline, $cpercent;
      } else {
        push @percentline, 0;
      }
    }
    
    my $percentline = join "\t", @percentline;
    
    push @percentout, $percentline;
  
  }
  
  $samplecount +=1;
  
  my $left = $samplecount - 0.2;
  
  my $right = $samplecount + 0.2;
  
  my $cumulpercent = 0;
  
  foreach my $s (@origcolors) {
    if (defined $localdip{$s}) {
      my $percent = sprintf "%.4f", (100*$localdip{$s}/$windowcount);
      if ($percent > 0) {
        my $start = $cumulpercent;
        $cumulpercent += $percent;
        print "rect($left,$start,$right,$cumulpercent,col=\"$s\",border=NA)\n";
      }
    }
  }
  
  my $result = join "\n", @out;
  unless ( open(OUT, ">$outfile") ) {
      print "Cannot open file \"$outfile\" to write to!!\n\n";
      exit;
  }
  print OUT $result;
  close (OUT);
  
  my $dipresult = join "\n", @outdip;
  unless ( open(DOUT, ">$dipoutfile") ) {
      print "Cannot open file \"$dipoutfile\" to write to!!\n\n";
      exit;
  }
  print DOUT $dipresult;
  close (DOUT);

}

print "\n\n\n";

foreach my $s (keys %fullsisters) {
  print "$s\t$fullsisters{$s}\n";
}

my $presult = join "\n", @percentout;
unless ( open(POUT, ">$percentoutfile") ) {
    print "Cannot open file \"$percentoutfile\" to write to!!\n\n";
    exit;
}
print POUT $presult;
close (POUT);
