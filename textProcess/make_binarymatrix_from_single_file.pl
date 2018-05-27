#!/usr/bin/perl
%rowlist=();
%collist=();
%flag=();
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	@t = split(/\t/);
	$id = $t[0].":".$t[1];
	$flag{$id} = 1;
	$rowlist{$t[0]}=1;
	$collist{$t[1]}=1;
	}
close(IN);

print "ID";
foreach $c(sort keys %collist){
	print "\t$c";
	}
print "\n";

foreach $r(keys %rowlist){
	print "$r";
	foreach $c(sort keys %collist){
		$id = $r.":".$c;
		if(not exists $flag{$id}){$flag{$id} = 0;}
		print "\t$flag{$id}";
		}
	print "\n";
	}
