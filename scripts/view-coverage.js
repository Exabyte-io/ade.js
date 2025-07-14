#!/usr/bin/env node

const { exec } = require("child_process");
const path = require("path");
const fs = require("fs");

const coveragePath = path.join(__dirname, "..", "coverage", "index.html");

if (!fs.existsSync(coveragePath)) {
    console.error('❌ Coverage report not found. Run "npm run test:coverage:html" first.');
    process.exit(1);
}

console.log("📊 Opening coverage report in browser...");

// Open coverage report in default browser
let command;
if (process.platform === "win32") {
    command = "start";
} else if (process.platform === "darwin") {
    command = "open";
} else {
    command = "xdg-open";
}

exec(`${command} "${coveragePath}"`, (error) => {
    if (error) {
        console.error("❌ Failed to open coverage report:", error.message);
        console.log(`📁 Coverage report is available at: ${coveragePath}`);
        process.exit(1);
    }
    console.log("✅ Coverage report opened successfully!");
});
