## Structure

Simulations of the Brownian bug model are run in C++ in `simulation` and output files are then treated separately with R to produce figures, in `figure`.

## Reference platform

Results presented in the article were obtained with the following configuration:

**Platfom**
```
Linux #137~16.04.1-Ubuntu
```

**Hardware information**

```
flags        : fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc art arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc cpuid aperfmperf pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand lahf_lm abm 3dnowprefetch cpuid_fault epb invpcid_single pti ssbd ibrs ibpb stibp tpr_shadow vnmi flexpriority ept vpid fsgsbase tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm mpx rdseed adx smap clflushopt intel_pt xsaveopt xsavec xgetbv1 xsaves dtherm ida arat pln pts hwp hwp_notify hwp_act_window hwp_epp md_clear flush_l1d
model name   : Intel(R) Xeon(R) CPU E3-1240 v5 @ 3.50GHz
vendor_id    : GenuineIntel
```

**Compilation information for C++ and R**

Simulation were run with C++, with the following configuration:

```
g++ (Ubuntu 5.4.0-6ubuntu1~16.04.12) 5.4.0 20160609
```

with compilation flags -O3 -std=c++11 and library GSL v.2.6.
