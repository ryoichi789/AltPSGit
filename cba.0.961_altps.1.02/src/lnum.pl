#!/usr/bin/perl

for (my $n = 1; <>; ++$n) {
	s/\r//g; # for Windows->Unix
	s/\t/    /g; # tab -> 4 spaces
	s/^\/\/\s*\S+:\s*\n$/\/\/\n/;
	printf(STDOUT "%4d: %s", $n, $_);
}

1;
