/**
 * mock-fastq-generator - UI Interactions & App Logic
 */

document.addEventListener("DOMContentLoaded", () => {
    let phredChart = null;

    // Elements
    const form = document.getElementById("configForm");
    const previewBtn = document.getElementById("previewBtn");
    const toggleScaleBtn = document.getElementById("toggleScaleBtn");
    const generateBtn = document.getElementById("generateBtn");
    const statusMsg = document.getElementById("statusMessage");
    const scoreInputMode = document.getElementById("scoreInputMode");
    const minValInput = document.getElementById("minVal");
    const maxValInput = document.getElementById("maxVal");
    const minValLabel = document.getElementById("minValLabel");
    const maxValLabel = document.getElementById("maxValLabel");

    let isAbsolutePhred = true;

    // Initialize chart with empty data
    initChart();
    
    // Draw initial preview based on default values
    updatePreview();

    // Event Listeners
    previewBtn.addEventListener("click", () => {
        updatePreview();
    });

    toggleScaleBtn.addEventListener("click", () => {
        isAbsolutePhred = !isAbsolutePhred;
        toggleScaleBtn.textContent = isAbsolutePhred ? "Show ASCII Scale" : "Show Absolute Phred";
        updatePreview();
    });

    form.addEventListener("keydown", (e) => {
        if (e.key === "Enter") {
            e.preventDefault();
            updatePreview();
        }
    });

    scoreInputMode.addEventListener("change", () => {
        const isPhred = scoreInputMode.value === "phred";
        
        // Update labels
        minValLabel.innerHTML = isPhred 
            ? 'Min Score <span class="info-icon" title="The lowest allowed quality score.">ℹ️</span>'
            : 'Min ASCII Code <span class="info-icon" title="The lowest allowed quality score, encoded as ASCII. 40 equals Phred 7.">ℹ️</span>';
            
        maxValLabel.innerHTML = isPhred
            ? 'Max Score <span class="info-icon" title="The highest allowed quality score at the peak.">ℹ️</span>'
            : 'Max ASCII Code <span class="info-icon" title="The highest allowed quality score at the peak. 73 equals Phred 40.">ℹ️</span>';

        // Convert the current values in the input boxes
        let currentMin = parseInt(minValInput.value, 10);
        let currentMax = parseInt(maxValInput.value, 10);
        
        if (!isNaN(currentMin)) {
            minValInput.value = isPhred ? currentMin - 33 : currentMin + 33;
        }
        if (!isNaN(currentMax)) {
            maxValInput.value = isPhred ? currentMax - 33 : currentMax + 33;
        }
        
        updatePreview();
    });

    generateBtn.addEventListener("click", () => {
        statusMsg.textContent = "Generating dataset...";
        generateBtn.disabled = true;

        // Use setTimeout to allow UI to update before heavy synchronous JS computation
        setTimeout(() => {
            try {
                generateAndDownload();
                statusMsg.textContent = "Success! Files downloaded.";
                statusMsg.style.color = "#4ade80"; // Green
            } catch (err) {
                console.error(err);
                statusMsg.textContent = "Error: " + err.message;
                statusMsg.style.color = "#f87171"; // Red
            } finally {
                generateBtn.disabled = false;
            }
        }, 50);
    });

    // Live updates on form input changes
    form.addEventListener("input", () => {
        updatePreview();
    });

    const decayModelSelect = document.getElementById("decayModel");
    const decayRateGroup = document.getElementById("decayRateGroup");
    const decayRateInput = document.getElementById("decayRate");
    const decayRateLabel = document.getElementById("decayRateLabel");
    const centerGroup = document.getElementById("centerGroup");
    const centerInput = document.getElementById("center");
    const centerLabel = document.getElementById("centerLabel");
    const stdDevGroup = document.getElementById("stdDevGroup");

    decayModelSelect.addEventListener("change", () => {
        const model = decayModelSelect.value;
        if (model === "exponential") {
            decayRateGroup.style.display = "block";
            stdDevGroup.style.display = "none";
            centerGroup.style.display = "none"; // Exponential model does not use center (C)
            decayRateLabel.innerHTML = 'Decay Rate (λ) <span class="info-icon" title="Exponential decay rate. Sub-exponential initial plateau via x^1.25.">ℹ️</span>';
            if (decayRateInput.value === "1" || decayRateInput.value === "1.0") {
                decayRateInput.value = "0.0008";
            }
        } else if (model === "sigmoidal") {
            decayRateGroup.style.display = "block";
            stdDevGroup.style.display = "none";
            centerGroup.style.display = "block";
            decayRateLabel.innerHTML = 'Drop Steepness (k) <span class="info-icon" title="Sigmoidal drop steepness. Higher values create sharp step-function drops (e.g. k=1.0).">ℹ️</span>';
            centerLabel.innerHTML = 'Drop Position (C) <span class="info-icon" title="Nucleotide coordinate where quality drops to the midpoint.">ℹ️</span>';
            if (decayRateInput.value === "0.0008") {
                decayRateInput.value = "1.0";
            }
        } else {
            // gaussian
            decayRateGroup.style.display = "none";
            stdDevGroup.style.display = "block";
            centerGroup.style.display = "block";
            centerLabel.innerHTML = 'Peak Center (C) <span class="info-icon" title="Base position where read quality is highest.">ℹ️</span>';
        }
        updatePreview();
    });

    // Trigger initially
    decayModelSelect.dispatchEvent(new Event("change"));

    function getParamsFromForm() {
        const mode = document.getElementById("scoreInputMode").value;
        let min = parseInt(document.getElementById("minVal").value, 10);
        let max = parseInt(document.getElementById("maxVal").value, 10);
        
        // Convert to ASCII if user input was Absolute Phred
        if (mode === "phred") {
            min += 33;
            max += 33;
        }

        return {
            templateSequence: document.getElementById("templateSequence").value,
            upstreamSequence: document.getElementById("upstreamSequence").value,
            leftMargin: parseInt(document.getElementById("leftMargin").value, 10),
            totalLength: parseInt(document.getElementById("totalLength").value, 10),
            numberOfSequences: parseInt(document.getElementById("numberOfSequences").value, 10),
            center: parseInt(document.getElementById("center").value, 10),
            minVal: min,
            maxVal: max,
            stdDev: parseInt(document.getElementById("stdDev").value, 10),
            noiseLevel: parseFloat(document.getElementById("noiseLevel").value),
            pairedEnd: document.getElementById("pairedEnd").checked,
            decayModel: document.getElementById("decayModel").value,
            decayRate: parseFloat(document.getElementById("decayRate").value),
            binnedQuality: document.getElementById("binnedQuality").checked,
            homopolymerPenalty: document.getElementById("homopolymerPenalty").checked
        };
    }

    function updatePreview() {
        const params = getParamsFromForm();
        // Just generate 1 sequence to get the plot data
        params.numberOfSequences = 1; 
        const result = assembleSequences(params);
        
        if (result.plotData) {
            renderChart(result.plotData);
        }
    }

    function generateAndDownload() {
        const params = getParamsFromForm();
        const result = assembleSequences(params);
        
        // Download R1
        const r1String = recordsToFastqString(result.r1Records);
        downloadFile("dataset_R1.fastq", r1String);
        
        // Download R2 if paired
        if (params.pairedEnd) {
            const r2String = recordsToFastqString(result.r2Records);
            downloadFile("dataset_R2.fastq", r2String);
        }

        // Update plot with the first sequence from this actual run
        if (result.plotData) {
            renderChart(result.plotData);
        }
    }

    function downloadFile(filename, content) {
        const blob = new Blob([content], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }

    function initChart() {
        const ctx = document.getElementById('phredChart').getContext('2d');
        
        // Gradient fill
        const gradient = ctx.createLinearGradient(0, 0, 0, 400);
        gradient.addColorStop(0, 'rgba(59, 130, 246, 0.5)'); // Blue
        gradient.addColorStop(1, 'rgba(59, 130, 246, 0.0)');
        
        phredChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: [],
                datasets: [{
                    label: 'Phred Score (ASCII Value)',
                    data: [],
                    borderColor: '#3b82f6',
                    backgroundColor: gradient,
                    borderWidth: 2,
                    pointRadius: 0,
                    pointHoverRadius: 4,
                    fill: true,
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: {
                        labels: { color: '#f8fafc' }
                    }
                },
                scales: {
                    x: {
                        grid: { color: 'rgba(255, 255, 255, 0.1)' },
                        ticks: { color: '#94a3b8' },
                        title: {
                            display: true,
                            text: 'Position (Index)',
                            color: '#94a3b8'
                        }
                    },
                    y: {
                        grid: { color: 'rgba(255, 255, 255, 0.1)' },
                        ticks: { color: '#94a3b8' },
                        title: {
                            display: true,
                            text: 'ASCII Score',
                            color: '#94a3b8'
                        },
                        min: 0,
                        max: 126
                    }
                }
            }
        });
    }

    function renderChart(numericScores) {
        // Transform data based on selected scale
        const data = isAbsolutePhred ? numericScores.map(score => score - 33) : numericScores;
        const labels = Array.from({length: data.length}, (_, i) => i);
        
        phredChart.data.labels = labels;
        phredChart.data.datasets[0].data = data;
        
        // Update Chart Labels
        phredChart.data.datasets[0].label = isAbsolutePhred ? 'Absolute Phred Score' : 'Phred Score (ASCII Value)';
        phredChart.options.scales.y.title.text = isAbsolutePhred ? 'Phred Score' : 'ASCII Score';
        
        // Adjust Y axis min/max dynamically to make graph look better
        const minScore = Math.min(...data);
        const maxScore = Math.max(...data);
        
        // Set min to minScore - 10, but don't let Absolute Phred go below 0
        phredChart.options.scales.y.min = isAbsolutePhred ? Math.max(0, minScore - 10) : Math.max(33, minScore - 10);
        phredChart.options.scales.y.max = maxScore + 10;
        
        phredChart.update();
    }
});
