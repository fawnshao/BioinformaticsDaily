#!/usr/bin/perl
if(@ARGV<4){die "Usage: perl $0 <file1> <file2> <order of the key in file1> <order of the key in file2> > <output>\nthe order should be 0-based\n\n";}
%hash=();

open(IN,$ARGV[0]);
print STDERR "reading $ARGV[0]\n";
while(<IN>){
        chomp;
        @t=split(/\t/);
        $id=$t[$ARGV[2]];
        $hash{$id}=$hash{$id}.$_."\n";
        }
close(IN);
print STDERR "finished $ARGV[0] hash\n";

$size = scalar(@t) - 1;
$na = "/\t"x$size."/";

open(IN,$ARGV[1]);
print STDERR "reading $ARGV[1]\n";
while(<IN>){
        chomp;
        @t=split(/\t/);
        if(not exists $hash{$t[$ARGV[3]]}){print "$_\t$na\n";}
        else{
                @tt=split(/\n/,$hash{$t[$ARGV[3]]});
                foreach $line(@tt){
                        print "$_\t$line\n";
                        }
                }
        }
close(IN);
print STDERR "finished $ARGV[1] hash\n";