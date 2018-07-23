#!/usr/bin/perl
%hash = ();
open(IN, $ARGV[0]) or die "can not open $ARGV[0]\n";
while(<IN>){
	chomp;
	@t = split(/\t/);
	$hash{$t[0]} .= $t[1].";";
}
close(IN);
foreach $id(sort keys %hash){
	my @t = split(/;/, $hash{$id});
	my @ut = keys { map { $_ => 1 } @t };
	my $a = join(";", @ut);
	print "$id\t$a\n";
}