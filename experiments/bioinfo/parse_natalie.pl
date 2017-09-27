#!/usr/bin/perl

open F1, "Mus_musculus.xml";
open F2, "Homo_sapiens.xml";

my %map1;
my %map2;

my $t;

$t = 0;
while(<F1>) {
    chomp;
    if (/data key.*>(.*)</) {
	$map1{$1} = $t++;
    }
}

$t = 0;
while(<F2>) {
    chomp;
    if (/data key.*>(.*)</) {
        $map2{$1} = $t++;
    }
}

open F3, "Mus_musculus_Homo_sapiens_ortholog_candidates_0.4.txt";
<F3>;

my $i;
my $j;

while (<F3>) {
    chomp;
    my @ids = split;
#    next unless defined $map1{$ids[0]};
    $i = $map1{$ids[0]};
    for ($k = 1; $k <= $#ids; $k++) {
	$j = $map2{$ids[$k]};
	print "$i $j\n";
    }
}
