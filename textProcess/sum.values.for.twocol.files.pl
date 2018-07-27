#!/usr/bin/perl
### input files has two columns, key - values;
%sum = ();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	$sum{$t[0]} += $t[1];
}
close(IN);

foreach $gene(sort keys %sum){
	print "$gene\t$sum{$gene}\n";
}