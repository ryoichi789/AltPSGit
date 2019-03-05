#!/usr/bin/perl

for (<>) {
	s/\r//g; # for Windows->Unix
	print $_;
}

1;
