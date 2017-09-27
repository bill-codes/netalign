#!/usr/bin/env PERL_UNICODE=S perl

open (BIBLEEN, "<utf8", "bib.EN");
open (BIBLEES, "<utf8", "bib.ES");
open (WORDEN, "<utf8", "word.EN");
open (WORDES, "<utf8", "word.ES");
open (OUTFILE, ">L.mat");

my $op;
my $verseName;
my $verseText;
my @vn1;
my @vn2;
my @list1;
my @list2;
my %worden;
my %wordes;
my $n1 = 0;
my $n2 = 0;

while (<WORDEN>) {
    chomp;
    $worden{$_} = $n1++;
}

$n = 0;

while (<WORDES>) {
    chomp;
    $wordes{$_} = $n2++;
}

sub parser {
    ($text, $hash) = @_;
    $text =~ s/[?;:!,.'"]//g;
    $text =~ tr/[A-Z]/[a-z]/;
    my @words = split /\s+/, $text;
    my @list;
    foreach $word (@words) {
	push (@list, $hash->{$word}) if exists $hash->{$word};
    }
    return \@list;
}

while (<BIBLEEN>) {
    chomp;
    if ($op == 0) {
	if (/<seg id="(.*?)".*?>(.*)/) {
	    $op = 1;
	    $verseName = $1;
	    $verseText = $2;
	    my $listpointer = &parser($verseText, \%worden);
	    push @vn1, $verseName;
	    push @list1, $listpointer;
	}
    }
    else {
	if (/<\/seg>/) {
	    $op = 0;
	}
	else {
	    $verseText = $verseText.$_;
	}
    }
}


while (<BIBLEES>) {
    chomp;
    if ($op == 0) {
	if (/<seg id="(.*?)".*?>(.*)/) {
	    $op = 1;
	    $verseName = $1;
	    $verseText = $2;
	    my $listpointer = &parser($verseText, \%wordes);
	    push @vn2, $verseName;
	    push @list2, $listpointer;
	}
    }
    else {
	if (/<\/seg>/) {
	    $op = 0;
	}
	else {
	    $verseText = $verseText.$_;
	}
    }
}

my $edges = {};
my $nedges;

for $t (0..$n1-1) {
    $edges->{$t} = {};
}

for $t (0..$#vn1) {
    foreach $i (@{$list1[$t]}) {
        $tmp = $edges->{$i};
        foreach $j (@{$list2[$t]}) {
            if (not exists $tmp->{$j}) {
                $tmp->{$j} = 1;
#                $nedges++;
            }
            else {
                $nedges++ if $tmp->{$j} == 1;
                $tmp->{$j}++;
            }
        }
    }
}

print OUTFILE "$n1 $n2 $nedges\n";
for $t (0..$n1-1) {
    foreach (keys %{$edges->{$t}}) {
        if ($edges->{$t}->{$_} > 1) {
            if ($edges->{$t}->{$_} <= 25) {
                my $tmp = $edges->{$t}->{$_} / 25;
                print OUTFILE "$t $_ $tmp\n";
            }
            else {
                print OUTFILE "$t $_ 1\n";
            }
        }
    }
}
