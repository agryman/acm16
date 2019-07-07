2019-07-06

G'day Boss,

I am indeed installed in Melbourne: I was thinking about dropping
in on you on my way here, but I wasn't given much notice after
we'd finally got the visa sorted out - and my boss wanted me
here quickly so that I could do some of his teaching!


Hello Arthur,

it's good to welcome you to the project and to celebrate
your escape from the dark side - although maybe I am
now in the twilight myself, being back doing pure maths.

Anyway, I hope my ACM webpages have answered some of your
questions. There are indeed various files of pre-calculated
(SO(3)-reduced) SO(5) CG coefficients. You can download
these files from the webpages. Calculating them on the fly
would be too time-consuming. The 2009 paper with Caprio
shows how they can be calculated. However, I used a
C program which I wrote before that: it isn't particularly
efficient, but I could happily let it calculate in the
background while I worked on something else.

Some information on how the ACM code accesses these 
coefficients is given in section 6.1 of the 2016 paper 
I wrote with David
(Computer Physics Communications 200 (2016), 220-253).
In the ACM code, the procedure
      load_CG_table()
is used to read one data file and install its values
into memory (into a Maple table named CG_coeffs).

As for the main ACM code itself, have you obtained the 
published version (1.4)? Earlier versions certainly had 
some bugs, weren't so well tested, and weren't so well 
documented. I remember that my efforts to eke out one 
bug were held back by my use of the "remember" option 
for some of the procedures. It might be a good idea to
take out all the "remember"s, if you're testing a new
language version against the 1.4 version.

I wonder which language you're intending to translate
it into. You'll of course need access to some very
efficient linear algebra routines that enable large
floating point matrices (of size a few hundred)
to be diagonalised. This translation could be a lot
of work - I hope that it will be worthwhile!

Cheers,
Trevor

