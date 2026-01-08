# Installing Mosek for R (Step-by-Step Manual)

## 0. Prerequisites checklist

- Windows 11 64-bit

- R 64-bit (check in R: R.version$arch should be "x86_64")

- A working MOSEK installation (Optimizer Suite)

- Rtools matching your R version (needed because the interface is built/installed via MOSEK’s builder on Windows) 

## 1. Install MOSEK Optimization Suite (MSI)

Download the Windows 64-bit x86 MOSEK installer (MSI) from MOSEK downloads and run it. 

After installation, ensure the MOSEK binaries folder is on your PATH. MOSEK’s install guide explicitly calls this out as necessary so Windows can find MOSEK shared libraries (DLLs). 

The folder is typically like:

```
<MSKHOME>\mosek\11.0\tools\platform\win64x86\bin
```

where <MSKHOME> is your MOSEK install root. 

Quick test (Command Prompt / PowerShell)

Run: ```mosek```

You should see the MOSEK version banner (like you did earlier). If Windows says it can’t find mosek, your PATH isn’t set correctly. 

## 2. Install / activate the MOSEK license

Make sure MOSEK itself is licensed (academic license is fine). You can place the license file in MOSEK’s expected location or configure the environment variable MOSEK uses for the license. (If mosek runs and MOSEK can solve examples, you’re usually fine.)

Academic lincense can be requested at
https://www.mosek.com/products/academic-licenses/
You will receive an email license file and with further instructions.

If you later get “license” errors in Rmosek, it means MOSEK runs but the license location isn’t visible from that R session.


## 3. Install Rtools (Windows)

Install the Rtools version matching your installed R (e.g., R 4.5 → Rtools45). MOSEK’s builder has helpers to resolve Rtools on Windows, but Rtools still needs to be installed. 

Verify Rtools is visible

In R:

```
Sys.which("make")
Sys.which("gcc")
```

These should return non-empty paths.

## 4. Install the Rmosek interface using MOSEK’s builder.R

### 4.1 Locate builder.R

Find your MOSEK “Rmosek directory” (<RMOSEKDIR>). In MOSEK 11 it is typically inside your MOSEK install tree, e.g. something like:

```
<MSKHOME>\mosek\11.0\tools\platform\win64x86\rmosek
```

(Exact path can differ; the key is: the directory that contains builder.R.) 


### 4.2 Run builder from inside R

In R:

```
source("<RMOSEKDIR>/builder.R")
attachbuilder(what_mosek_bindir = "<MSKHOME>/MOSEK/11.0/tools/platform/win64x86/bin")
install.rmosek()
```

## 5. Verify from R 

In a fresh R session:

```
library(Rmosek)
Rmosek::mosek_version()
```

Then try a tiny solve (MOSEK docs include examples; any small LP/SOCP test is fine).





