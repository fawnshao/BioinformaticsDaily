#!/usr/bin/perl
open(IN,$ARGV[0]);
%hash = ();
%strand = ();
while(<IN>){
	chomp;
	@t = split(/\t/);
	$id = $t[0]."\t".$t[1]."\t".$t[2];
	$strand{$id} = $t[5];
	if(not exists $hash{$id}){
		$hash{$id} = $t[3];
		}
		else{
			$hash{$id} = $hash{$id}.";".$t[3];
			}
	}
close(IN);

foreach $a(sort keys %hash){
	print "$a\t$hash{$a}\t.\t$strand{$a}\n";
	}
