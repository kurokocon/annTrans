#!/usr/bin/env perl

use strict;
use warnings;
my $i = 0;
while (<>) {
	if (/>(\S+)/) {
		my $acc = $1;
		print "NA\t$acc\n";
		$i++;
	}
}
