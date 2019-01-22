#!/usr/bin/perl
%hash = ();
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	if(not exists $hash{$t[0]} or $hash{$t[0]} < $t[1]){
		$hash{$t[0]} = $t[1];
	}
}
close(IN);
foreach $id(keys %hash){
	print "$id\t$hash{$id}\n";
}