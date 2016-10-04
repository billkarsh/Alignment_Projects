
## Bills_aligner_tools.py

**Contributed by Joergen Kornfeld with these comments:**

```
I added all python stuff that I used for the analysis and the final stack rendering to github and
if you want, you could pull it into your main repository. The tile rendering is pretty fast IMO
(11TB in <24h on 32 CPU cores, basically IO limited), uses qsub and has many options that were
useful for me (different image normalization methods, gaussian filtering for very noisy images,
gray-value inversion, border cropping etc).

However, it is not really end-user polished and easy to use at this point, so if you know anybody
that might be interested in using it, just direct them to me.

Cheers,
Jorgen

--
Joergen Kornfeld
Max-Planck-Institute for Medical Research Jahnstrasse 29
69120 Heidelberg, Germany
```

_fin_
