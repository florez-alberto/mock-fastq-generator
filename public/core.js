/**
 * mock-fastq-generator - Core Logic (JavaScript Port)
 * This replicates the logic found in src/mock_fastq_generator/core.py
 */

const ASCII_BASE = 33;
const BASES = ['A', 'C', 'G', 'T'];

// Helper for normally distributed random numbers (Box-Muller transform)
function randomNormal(mean = 0, stdDev = 1) {
    let u = 1 - Math.random(); // Converting [0,1) to (0,1]
    let v = Math.random();
    let z = Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
    return z * stdDev + mean;
}

// Generate an Illumina-style header
function generateReadHeader() {
    const randomChars = () => Math.random().toString(36).substring(2, 6).toUpperCase();
    const randomDigits = () => Math.floor(1000 + Math.random() * 9000).toString();
    return `@${randomChars()}:${randomDigits()}`;
}

// Reverse complement a DNA sequence
function reverseComplement(sequence) {
    const complementMap = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N' };
    return sequence
        .toUpperCase()
        .split('')
        .reverse()
        .map(base => complementMap[base] || base)
        .join('');
}

// Add substitution noise
function addNoiseToSequence(sequence, noiseLevel) {
    let noisy = '';
    for (let i = 0; i < sequence.length; i++) {
        if (Math.random() < noiseLevel) {
            noisy += BASES[Math.floor(Math.random() * BASES.length)];
        } else {
            noisy += sequence[i];
        }
    }
    return noisy;
}

// Generate Phred scores using the parametric models
function generatePhredScores(sequence, center, minVal, maxVal, stdDev, noiseLevel, decayModel = 'gaussian', decayRate = 0.0008, binnedQuality = false, homopolymerPenalty = false) {
    let scores = [];
    let asciiString = '';
    const numPoints = sequence.length;
    
    // Homopolymer calculation
    let homopolymerMask = new Array(numPoints).fill(false);
    if (homopolymerPenalty) {
        let runLen = 1;
        for (let i = 1; i < numPoints; i++) {
            if (sequence[i] === sequence[i-1]) {
                runLen++;
                if (runLen > 4) {
                    homopolymerMask[i] = true;
                }
            } else {
                runLen = 1;
            }
        }
    }

    for (let x = 0; x < numPoints; x++) {
        let curve = 0;
        
        if (decayModel === 'exponential') {
            curve = (maxVal - minVal) * Math.exp(-decayRate * Math.pow(x, 1.25)) + minVal;
        } else if (decayModel === 'sigmoidal') {
            curve = minVal + (maxVal - minVal) / (1 + Math.exp(decayRate * (x - center)));
        } else {
            // gaussian
            let exponent = -Math.pow((x - center), 2) / (2 * Math.pow(stdDev, 2));
            let gaussian = Math.exp(exponent);
            curve = gaussian * (maxVal - minVal) + minVal;
        }
        
        // Add noise directly proportional to noiseLevel parameter
        let noise = randomNormal(0, noiseLevel);
        let score = curve + noise;
        
        if (homopolymerMask[x]) {
            score -= 10;
        }
        
        // Clip to bounds
        score = Math.max(minVal, Math.min(maxVal, score));
        score = Math.round(score);
        
        if (binnedQuality) {
            // NovaSeq bins in Absolute Phred: 2, 12, 23, 37
            // minVal/maxVal are passed in as ASCII, so we add 33 to bins
            const bins = [2 + ASCII_BASE, 12 + ASCII_BASE, 23 + ASCII_BASE, 37 + ASCII_BASE];
            let closest = bins[0];
            let minDiff = Math.abs(score - bins[0]);
            for (let b = 1; b < bins.length; b++) {
                let diff = Math.abs(score - bins[b]);
                if (diff < minDiff) {
                    minDiff = diff;
                    closest = bins[b];
                }
            }
            score = closest;
        }
        
        scores.push(score);
        asciiString += String.fromCharCode(score);
    }
    
    return { asciiString, numericScores: scores };
}

// Main assembly function (returns FASTQ records)
function assembleSequences(params) {
    const {
        templateSequence,
        upstreamSequence,
        leftMargin,
        totalLength,
        numberOfSequences,
        center,
        minVal,
        maxVal,
        stdDev,
        noiseLevel,
        pairedEnd,
        decayModel = 'gaussian',
        decayRate = 0.0008,
        binnedQuality = false,
        homopolymerPenalty = false
    } = params;

    let r1Records = [];
    let r2Records = [];

    // Clean input sequences
    const template = templateSequence.replace(/\s+/g, '').toUpperCase();
    const upstream = upstreamSequence.replace(/\s+/g, '').toUpperCase();

    // Calculate margins
    const requiredLength = upstream.length + leftMargin + template.length;
    let rightMargin = 0;
    if (totalLength > requiredLength) {
        rightMargin = totalLength - requiredLength;
    }

    // We collect the numerical scores of the FIRST generated read to show in the plot
    let plotData = null;

    for (let i = 0; i < numberOfSequences; i++) {
        // Build sequence parts
        let leftSeq = Array.from({ length: leftMargin }, () => BASES[Math.floor(Math.random() * BASES.length)]).join('');
        let rightSeq = Array.from({ length: rightMargin }, () => BASES[Math.floor(Math.random() * BASES.length)]).join('');
        
        let noisyTemplate = addNoiseToSequence(template, noiseLevel);
        let finalSeq = upstream + leftSeq + noisyTemplate + rightSeq;
        
        // Ensure exact totalLength if requested
        if (finalSeq.length > totalLength) {
            finalSeq = finalSeq.substring(0, totalLength);
        } else if (finalSeq.length < totalLength) {
            let padding = Array.from({ length: totalLength - finalSeq.length }, () => 'N').join('');
            finalSeq += padding;
        }

        const seqLen = finalSeq.length;
        const header = generateReadHeader();
        const phredData = generatePhredScores(
            finalSeq, center, minVal, maxVal, stdDev, noiseLevel,
            decayModel, decayRate, binnedQuality, homopolymerPenalty
        );
        
        if (i === 0) {
            plotData = phredData.numericScores;
        }

        r1Records.push([
            pairedEnd ? `${header}/1` : header,
            finalSeq,
            "+",
            phredData.asciiString
        ]);

        if (pairedEnd) {
            const r2Seq = reverseComplement(finalSeq);
            const r2PhredData = generatePhredScores(
                r2Seq, center, minVal, maxVal, stdDev, noiseLevel,
                decayModel, decayRate, binnedQuality, homopolymerPenalty
            );
            r2Records.push([
                `${header}/2`,
                r2Seq,
                "+",
                r2PhredData.asciiString
            ]);
        }
    }

    return { r1Records, r2Records, plotData };
}

// Convert record array to FASTQ formatted string
function recordsToFastqString(records) {
    return records.map(record => record.join('\n')).join('\n') + '\n';
}
