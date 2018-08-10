# calculate the effective length for each anchor from bedtools intersect files
#!/usr/bin/perl
%hash = ();
open(IN, $ARGV[0]);
while(<IN>){
    chomp;
    @t = split(/\t/);
    if($t[11] < 1000) {
    	$hash{$t[3]} += $t[10];
    }
    else{
    	$hash{$t[3]} += 1000;
    }
}
close(IN);

foreach $anchor(keys %hash){
	print "$anchor\t$hash{$anchor}\n";
}
