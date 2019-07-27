#!/usr/bin/perl
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
%chrs = ();
%hleft = ();
%hright = ();
%strands = ();
while(<IN>){
	chomp;
	@t = split(/\t/);
	if(not exists $htss{$t[3]}){
		$chrs{$t[3]} = $t[0];
		$hleft{$t[3]} = $t[1];
		$hright{$t[3]} = $t[2];
		$strands{$t[3]} = $t[5];
	}
	else{
		if($t[1] > $hleft{$t[3]}){
			$hleft{$t[3]} = $t[1];
		}
		if($t[2] < $hright{$t[3]}){
			$hright{$t[3]} = $t[2];
		}
	}
}
close(IN);

foreach $id (keys %chrs){
	print "$chrs{$id}\t$hleft{$id}\t$hright{$id}\t$id\t0\t$strands{$id}\n";
}