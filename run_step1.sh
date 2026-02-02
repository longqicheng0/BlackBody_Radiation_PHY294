#!/bin/bash
# Quick setup and run script for Step 1 analysis

echo "======================================"
echo "Blackbody Lab Step 1 Analysis Setup"
echo "======================================"
echo ""

# Check Python version
if command -v python3 &> /dev/null; then
    PYTHON_CMD=python3
    echo "✓ Found python3"
elif command -v python &> /dev/null; then
    PYTHON_CMD=python
    echo "✓ Found python"
else
    echo "✗ Python not found! Please install Python 3.9+"
    exit 1
fi

# Show Python version
echo "  Python version: $($PYTHON_CMD --version)"
echo ""

# Install dependencies if needed
echo "Installing dependencies..."
$PYTHON_CMD -m pip install --quiet --upgrade pip
$PYTHON_CMD -m pip install --quiet -r requirements.txt

if [ $? -eq 0 ]; then
    echo "✓ Dependencies installed"
else
    echo "✗ Failed to install dependencies"
    exit 1
fi

echo ""
echo "======================================"
echo "Running Step 1 Analysis..."
echo "======================================"
echo ""

# Run the analysis
$PYTHON_CMD -m bb_step1 "$@"
