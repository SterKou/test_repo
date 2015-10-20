#! /usr/bin/perl -w

use strict;

my %result    = ();
my $LIMIT     = 230;
my $AMP_LIMIT = 150;
my $TM_DOWN   = 58;
my $TM_UP     = 62;

my $rev = 0;

while(<>) { 

   chomp;
   s/\r//;
   my @row = split/\t/;

   my $id      = $row[0];
   my $seq     = $row[2];
   my $motif   = $row[1];  
   my $mot_len = length($motif); 

   #print STDERR "$id\n$seq\n$motif\n";

   if ($seq !~ /$motif/) {
      $motif = rev_com($motif);    
      #print STDERR "REV_COM: $motif\n";
      $rev++;

      if ($seq !~ /$motif/) {
         print STDERR "There is no $seq in $id\n";
         next;
      }
   }
   
   my $pos     = index($seq, $motif);
   my $start   = $pos-$LIMIT;
   my $end     = $pos+$mot_len;

   #print STDERR "Position: $pos\n";
 
   my $region_left   = substr($seq,$start,$LIMIT);   
   my $region_right  = substr($seq,$end,$LIMIT);
   $region_right     = reverse($region_right);
   $region_right     =~ tr/ATGCatgc/TACGtacg/;   

   #print STDERR "Start: $start\tEnd: $end\n";
   #print STDERR "$region_left\n$region_right\n";

   my $r_len    = length($region_left);
 
   for (my $i = 20; $i <=24; $i++) {  
      select_primer($id, $region_left, $r_len, "Forward", $i);
      select_primer($id, $region_right, $r_len, "Reverse", $i);
   }

}

print STDERR "Rev_complimented spacers: $rev\n";

for (keys %result) {

    my $id = $_;
    my %temp = %{$result{$id}};

   for (keys %temp) { 

     my $direct  = $_;
     my %temp_1  = %{$temp{$direct}};
     my $count   = 0;

      for ( sort{ $temp_1{$b}->[0] <=> $temp_1{$a}->[0] } keys %temp_1) {      
         $count++;       
         my $tm = $temp_1{$_}->[0];
         my $pos = $temp_1{$_}->[1];
    
         print "$id\t$direct\t$_\t$tm\t$pos\n" if($count <= 10);  
     }
     print STDERR "Total count for $id\t$count\n" if ($count < 10);

   }

}


sub rev_com {

   my $seq = shift;

  $seq = reverse($seq);
  $seq =~ tr/ATGCatgc/TACGtacg/;

  return $seq; 

}

sub select_primer {

   my ($id, $seq, $prod_len, $direct, $win) = @_;
   my $count  = 0;

   for (my $i=0; $prod_len >=$AMP_LIMIT; $i++) {

      my $primer = substr($seq, $i, $win);

      my $str  = $i+1;
      my $en   = $i+$win;
      $prod_len= $LIMIT-$en;

      #my $gc  = get_GC($primer);
      my $tm   = get_Tm($primer);

       if ( ($tm >=$TM_DOWN) && ($tm <=$TM_UP)) {

         $result{$id}->{$direct}->{$primer}->[0]= $tm;
         $result{$id}->{$direct}->{$primer}->[1]= $prod_len;
   
       }      

      $count++;

   }
   #print STDERR "Total primers: $count\n";

}

sub get_GC {

  my $seq = shift;
  my $len = length($seq);  
  my $count = 0;

  $count++ while($seq=~ /[GCgc]/g);    
  my $fraction = sprintf("%.0f", ($count/$len)*100);
  return $fraction;

}


sub get_Tm {

  my $seq = shift;
  my $tm  = 0;
  my $at  = 0;
  my $gc  = 0;

  $at++ while($seq=~ /[ATat]/g);
  $gc++ while($seq=~ /[GCgc]/g);
  
  $tm = $at*2 + $gc*4;
  return $tm;

}

