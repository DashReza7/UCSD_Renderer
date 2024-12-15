# Define the list of files to compress
$filesToCompress = @(
    "../ResultImages/cornellSimple.png",
    "../ResultImages/cornellNEE.png"
#     "../ResultImages/cornellRR.png",
#     "../ResultImages/dragon.png"
)

# Define the output ZIP file path
$outputZipFile = "../ResultImages/homework3a.zip"

# Filter out non-existent files
$availableFiles = @()
foreach ($file in $filesToCompress) {
    if (Test-Path -Path $file) {
        $availableFiles += $file
    }
}

# If there are available files, compress them
if ($availableFiles.Count -gt 0) {
    Compress-Archive -Path $availableFiles -DestinationPath $outputZipFile -Force
    Write-Host "Available files compressed successfully to $outputZipFile"
} else {
    Write-Host "No available files to compress."
}
