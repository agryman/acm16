2019-07-07

Hi Arthur (and David),

much has happened since version 16! I hope you're still not using version
16 David!!

When we first put our manuscript on the arxiv (arxiv.org/abs/1408.3824),
I aligned the manuscript and code version numbers, and decided to start
a new scheme at version 1.1 . The version that was eventually published
in CPC was version 1.4 . I attach the zip file containing that version;
this you get if you download the program code from the CPC library
(except I have removed the SO(5)>SO(3) coefficients because altogether
they're too big for an email - not even my unimelb email will allow it.
You can acquire these instead via my webpages).
There are a few other useful things in the zip file, as you'll see.

As written, the ACM code does quite a bit symbolic computation, keeping
things algebraic as long as possible. Then, all the matrix elements
get converted to floating point before diagonalisation is carried out.
If all you want to do is get the results after diagonalisation,
much of the symbolic computation is unnecessary, and the conversion to
floating point values could be done earlier.
I'll leave that up to you.

Bye for now,
Trevor.