#!/usr/bin/env PERL_UNICODE=S perl

open (FILE, "sense.EN");
open (OUTGRAPHFILE, ">graph.EN.smat");
open (WORDFILE, ">:utf8", "word.EN");

my $n;
my %words;
my %edgemap;
my @edges;
my $nedges;

while(<FILE>) {
    chomp;
    my @word = split(/ +/);
    if (substr($word[0], 0, 1) eq "S") {
	my @array;
	for $i (1 .. $#word) {
	    $words{$word[$i]} = $n++ unless exists $words{$word[$i]};
	    push @array, $words{$word[$i]};
	}
	for $i (0 .. $#array) {
	    for $j ($0 .. $#array) {
		next if $i == $j;
		my @edge = ($array[$i], $array[$j]);
		if (not exists $edgemap{$array[$i] * 10000 + $array[$j]}) {
		    $edgemap{$array[$i] * 10000 + $array[$j]} = 1;
		    $nedges++;
		    push @edges, [@edge];
		}
	    }
	}
    }
}

print OUTGRAPHFILE "$n $n $nedges\n";

foreach $edge (@edges) {
    print OUTGRAPHFILE "$edge->[0] $edge->[1] 1\n";
}

sub hashAscending {
    $words{$a} <=> $words{$b}
}

foreach (sort hashAscending keys %words) {
    print WORDFILE "$_\n";
}
