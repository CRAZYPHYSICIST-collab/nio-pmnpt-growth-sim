# Usage:  powershell -ExecutionPolicy Bypass -File run.ps1 -Mode all
#  Modes: test | exact | virtual | nonlinear | all
param([string]$Mode = "all", [string]$OutDir = "$PSScriptRoot\output")
$src = Get-Content -Raw -Encoding UTF8 "$PSScriptRoot\NiOSim.cs"
$cp = New-Object System.CodeDom.Compiler.CompilerParameters
[void]$cp.ReferencedAssemblies.Add("System.dll")
[void]$cp.ReferencedAssemblies.Add("System.Core.dll")
[void]$cp.ReferencedAssemblies.Add("System.Drawing.dll")
$cp.CompilerOptions = "/optimize+"
$cp.GenerateInMemory = $true
Add-Type -TypeDefinition $src -CompilerParameters $cp
[NiOSim.Program]::Main(@($Mode, $OutDir))
