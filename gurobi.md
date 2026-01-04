# Installing Gurobi for R (Step-by-Step Manual)

This document explains how to install Gurobi and its R interface so that it can be used by the `OptimalDesign` package.

The procedure consists of three independent steps:

- Installing the Gurobi optimizer

- Installing a Gurobi license

- Installing and testing the R interface

All three steps are required.

## 1. Installing Gurobi Optimizer
### 1.1 Download Gurobi

Go to the official Gurobi website:

https://www.gurobi.com/downloads/

Choose the installer corresponding to your operating system.

Download and run the installer.

### 1.2 Verify installation

After installation, open a terminal / command prompt.

Run: `gurobi_cl`

If Gurobi is installed correctly, you should see output similar to:

`Gurobi Optimizer version X.Y.Z`
`Copyright (c) Gurobi Optimization, LLC`

If the command is not found, restart your computer and try again.

## 2. Installing a Gurobi License

Gurobi requires a license to run. Free academic licenses are available.

### 2.1 Obtain an academic license

Go to:

https://www.gurobi.com/academia/academic-program-and-licenses/

Request an Academic Named-User License.

You will receive a license key of the form:

xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx

### 2.2 Install the license

In a terminal / command prompt, run:

`grbgetkey xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx`


This creates a license file gurobi.lic.

By default, the file is stored in:

Windows:
`C:\Users\<username>\gurobi.lic`

macOS / Linux:
`~/gurobi.lic`


## 3. Installing the Gurobi R Interface
### 3.1 Start R

Start a fresh R session (RStudio or plain R is fine).

### 3.2 Locate the Gurobi R directory

Gurobi installs its R interface inside the Gurobi installation directory.

Typical locations:

Windows

`C:/gurobi1300/win64/R`

macOS

`/Library/gurobi1300/macos_universal2/R`


Linux

`/opt/gurobi1300/linux64/R`

(Version number 1300 may differ.)

### 3.3 Install the R package

In R, run (adjust path if needed):

```
install.packages(
  "C:/gurobi1300/win64/R/gurobi_13.0-0.zip",
  repos = NULL
)
```

On macOS/Linux:

```
install.packages(
  "/Library/gurobi1300/macos_universal2/R/gurobi_13.0-0.tar.gz",
  repos = NULL,
  type = "source"
)
```

### 3.4 Load the package
`library(gurobi)`


If this succeeds without error, the R interface is installed.

## 4. Testing the Gurobi–R Installation

Run the following minimal example in R:

```
library(gurobi)

model <- list(
  A = matrix(c(1, 1), nrow = 1),
  obj = c(1, 1),
  rhs = 1,
  sense = "=",
  modelsense = "min",
  lb = c(0, 0)
)

result <- gurobi(model)
print(result$x)
```

Expected output:

`[1] 1 0`


(or an equivalent optimal solution)

If this works, Gurobi is correctly installed and usable from R.

## 5. Using Gurobi with OptimalDesign

Once Gurobi is installed and tested, run 

`library(OptimalDesign)`


The package will automatically detect Gurobi and use it for MISOCP problems.

## 6. Common Problems and Solutions
- Gurobi cannot find license

Verify that gurobi.lic exists

Check environment variable:

`Sys.getenv("GRB_LICENSE_FILE")`

- library(gurobi) fails

Make sure you installed Gurobi’s own R package, not a CRAN package

Restart R and try again

- Gurobi works in terminal but not in R

Restart RStudio

Ensure R and Gurobi are both 64-bit

## 7. Summary Checklist

Before reporting an installation problem, verify:

gurobi_cl works in terminal

License installed and valid

library(gurobi) works in R

Minimal R example solves successfully

