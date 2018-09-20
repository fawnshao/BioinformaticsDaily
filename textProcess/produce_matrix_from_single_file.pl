#!/usr/bin/perl
%tflist=();
%genelist=();
%flag=();
open(IN,$ARGV[0]);
while(<IN>){
	chomp;
	@t=split(/\t/);
	$id=$t[0].":".$t[1];
	# $flag{$id}=join("\t",$t[2],$t[3]);
	$flag{$id}=$t[2];
	$tflist{$t[0]}=1;
	$genelist{$t[1]}=1;
	}
close(IN);

@tfs=();
foreach $g(keys %tflist){
	push(@tfs,$g);
	}

print "Gene";
foreach $tf(@tfs){
	print "\t$tf";
	# print ".score\t$tf";
	# print ".pvalue";
	}
print "\n";

foreach $g(keys %genelist){
	print "$g";
	foreach $tf(@tfs){
		$id=$tf.":".$g;
		# if(not exists $flag{$id}){$flag{$id}="na\tna";}
		if(not exists $flag{$id}){$flag{$id}="";}
		print "\t$flag{$id}";
		}
	print "\n";
	}
