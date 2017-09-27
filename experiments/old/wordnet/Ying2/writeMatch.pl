#!/usr/bin/env PERL_UNICODE=S perl

open (EN, "<utf8", "word.EN");
open (ES, "<utf8", "word.ES");

my @wen;
my @wes;

while (<EN>) {
    chomp;
    push @wen, $_;
}

while (<ES>) {
    chomp;
    push @wes, $_;
}

open (MATCH, "out.txt");

while (<MATCH>) {
    chomp;
    my @temp = split / +/;
    print "$wen[$temp[0]-1] $wes[$temp[1]-1]\n";
}
