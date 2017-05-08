#!/usr/bin/python

def f1(a, b, c, d):
	print a, b, c, d
	return a + b + c + d

x = 5;
y = 6;
z = 7;

fm = lambda m: f1(m,x,y,z)
print fm(4)












x = 5;
y = 6;
z = 7;

fm = @(m) f1(m,x,y,z);
fm(4)

function [retval] = f1(a, b, c, d)
	retval = a + b + c + d;
