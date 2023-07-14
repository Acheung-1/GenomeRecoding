var messages = document.getElementById("messages");
var inputDNA = document.getElementById("inputDNA");
var button = document.getElementById("button");

clearInput.addEventListener("click", function(){
    // document.getElementById("inputDNA").innerHTML = "''";
    document.getElementById("a").innerHTML = "Protein Sequence: ";
    document.getElementById("b").innerHTML = "Input DNA - Number of CpGs: ";
    document.getElementById("c").innerHTML = "Input DNA - GC Content (%): " ;
    document.getElementById("d").innerHTML = "Optimized CpG Locations: " ;
    document.getElementById("e").innerHTML = "Average distance between each CpG in optimized: ";
    document.getElementById("f").innerHTML = "Optimized DNA - Number of CpGs: ";
    document.getElementById("g").innerHTML = "Optimized DNA - GC Content (%): " ;
    document.getElementById("h").innerHTML = "Optimized DNA - Sequence: " ;
});

convert.addEventListener("click",function(e){
    var valid = check(e, inputDNA.value.toUpperCase());
    if (valid === false) {
        alert("Not valid sequence")
        return;
    }


    var inputProteinSequence = convertToProtein(e, inputDNA.value.toUpperCase());
    document.getElementById("a").innerHTML = "Protein Sequence: " + inputProteinSequence;
    
    var inputNumCpG = CpG_DinucleoFrequency(e, inputDNA.value.toUpperCase());
    document.getElementById("b").innerHTML = "Input DNA - Number of CpGs: " + inputNumCpG.join(', ');
    
    var inputGCcont = gcContent(e, inputDNA.value.toUpperCase());
    document.getElementById("c").innerHTML = "Input DNA - GC Content (%): " + inputGCcont;
    
    // Finding which CpGs are REQUIRED and can't be changed
    var proteinZZ = Array.apply(null, Array(inputNumCpG.length)).map(function (x, i) { return 0; })
    var nucIndex = Array.apply(null, Array(inputNumCpG.length)).map(function (x, i) { return 0; })
    var need = []

    for (var i = 0; i < inputNumCpG.length; i++) {
        proteinZZ[i] =inputNumCpG[i] % 3;
        nucIndex[i] = inputNumCpG[i];
    }
    
    for (var i = 0; i < proteinZZ.length; i++) {
        if (proteinZZ[i] == 1){
            need.push(nucIndex[i]);
        }
    }

    // Converting protein sequence into maximum CpGs available
    var maxCpGDNA = convertToDNA(e, inputProteinSequence);
    
    var maxCpGLoc = CpG_DinucleoFrequency(e, maxCpGDNA);

    var optimizedCpGLoc = optimizeSpacing(e, need, maxCpGLoc);

    // This is the optimized sequence!
    var optimizedCpGDNA = optimizeSequence(e, inputProteinSequence, optimizedCpGLoc);

    document.getElementById("d").innerHTML = "Optimized CpG Locations: " + optimizedCpGLoc.join(', ') ;
    var optimizedCpGAveDist = averageDistance(e, optimizedCpGLoc);
    document.getElementById("e").innerHTML = "Average distance between each CpG in optimized: " + optimizedCpGAveDist;

    // cgatcgcctcgatcg need should be 1,10; orig: 1,5,10,14; max: 1,5,8,10,14
    document.getElementById("f").innerHTML = "Optimized DNA - Number of CpGs: " + optimizedCpGLoc.length;
    
    var optimizedGCcont = gcContent(e, optimizedCpGDNA);
    document.getElementById("g").innerHTML = "Optimized DNA - GC Content (%): " + optimizedGCcont ;
    
    document.getElementById("h").innerHTML = "Optimized DNA - Sequence: " + optimizedCpGDNA ;

    var p = CpG_DinucleoFrequency(e, optimizedCpGDNA)

    // To check to see optimized DNA seq does have optimized CpG Location
    // for (i = 0; i<optimizedCpGLoc.length; i++) {
    //     console.log(optimizedCpGLoc[i] === p[i]);
    // }
    // console.log(optimizedCpGLoc.length)
    // console.log(need)
    // console.log(optimizedCpGLoc)
});

var nucToAmino = {
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
    'TAT':'Y', 'TAC':'Y', 'TAA':'#', 'TAG':'#',
    'TGT':'C', 'TGC':'C', 'TGA':'#', 'TGG':'W',
                
    'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                
    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
                
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
};
  
  
var aminoToNuc_AT = {
    'F':'TTT', 
    'S':'TCA',
    'Y':'TAT',
    'C':'TGT',
    'W':'TGG',
    'L':'TTG',
    'P':'CCA',
    'H':'CAT',
    'Q':'CAG',
    'R':'AGG',
    'I':'ATA',
    'M':'ATG',
    'T':'ACA',
    'N':'AAT',
    'K':'AAG',
    'V':'GTG',
    'A':'GCA',
    'D':'GAT',
    'E':'GAG',
    'G':'GGG',
    '#':'TAG',
};
  
var aminoToNuc_CG = {
    'F':'TTC', 
    'S':'TCG',
    'Y':'TAC',
    'C':'TGC',
    'W':'TGG',
    'L':'CTC',
    'P':'CCG',
    'H':'CAC',
    'Q':'CAG',
    'R':'CGG',
    'I':'ATC',
    'M':'ATG',
    'T':'ACG',
    'N':'AAC',
    'K':'AAG',
    'V':'GTC',
    'A':'GCG',
    'D':'GAC',
    'E':'GAG',
    'G':'GGC',
    '#':'TAG',
};
  
var aminoToNuc_C = {
    'F':'TTC', 
    'S':'TCC',
    'Y':'TAC',
    'C':'TGC',
    'W':'TGG',
    'L':'CTC',
    'P':'CCC',
    'H':'CAC',
    'Q':'CAG',
    'R':'CGG',
    'I':'ATC',
    'M':'ATG',
    'T':'ACC',
    'N':'AAC',
    'K':'AAG',
    'V':'GTC',
    'A':'GCC',
    'D':'GAC',
    'E':'GAG',
    'G':'GGC',
    '#':'TAG',
};


function check(e, input) {
    e.preventDefault();
    valid = false
    for (var i = 0; i < input.length; i+=1) {
        if (input[i] == "A" || input[i] == "T" || input[i] == "C" || input[i] == "G") {
            valid = true;
        } else {
            valid = false;
            break;
        }
    }
    return valid
}

function convertToProtein(e, DNASeq) {
    /*
    Takes input as a string of DNA sequences
    Outputs the protein sequence as as string
    */
    e.preventDefault();
    var proteinSeq = "";
  
    for (var i = 0; i < DNASeq.length; i += 3) {
        if (i + 3 <= DNASeq.length) {
            proteinSeq += nucToAmino[DNASeq.slice(i, i + 3)];
        }
    }
  
    return proteinSeq;
    // document.getElementById("edited").innerHTML = nucToAmino[DNASeq];
}

function convertToDNA(e, proteinSeq) {
    /*
    Takes input as a string of protein sequence
    Outputs DNA sequence as a string, to maximize the number of CpGs, using aminoToNuc_CG dictionary
    */
    e.preventDefault();
    var DNASeq;
    DNASeq = "";
  
    for (var i = 0; i < proteinSeq.length; i += 1) {
        DNASeq += aminoToNuc_CG[proteinSeq[i]];
    }
  
    return DNASeq;
}
  
function gcContent(e, DNAseq) {
    /*
    Takes input as a string of DNA sequence
    Groups the numbers of Cs and Gs
    Groups the numbers of As and Ts
    Outputs total Gs and Cs / (Total nucleotides)
    */
    e.preventDefault();
    var at, gc;
    gc = 0;
    at = 0;
  
    for (var i = 0; i < DNAseq.length; i += 1) {
        if (DNAseq[i] === "C" || DNAseq[i] === "G") {
            gc ++;
        } else {
            at ++;
        }
    }
  
    return gc / (at + gc) * 100;
}
  
function CpG_DinucleoFrequency(e, DNAseq) {
    /*
    Takes input as a string of DNA sequence
    Counts the numbers of Cs followed by Gs
    Outputs an array corresponding to the location of Gs
    */
    e.preventDefault();
    var CpGlocation, counts, m;
    m = 0;
    counts = 0;
    CpGlocation = [];

    while (m < DNAseq.length) {
        if (DNAseq[m] === "C") {
            if (m < DNAseq.length - 1) {
                if (DNAseq[m + 1] === "G") {
                    counts += 1;
                    m += 1;
                    CpGlocation.push(m);
                }
            }
        }
        m += 1
    }
    return CpGlocation;
}
  
  
function optimizeSpacing(e, requiredCpG, CpGArray) {
    /*
    This function will has two inputs:
     necessaryCpG: An array containing the necessary locations required to be in the final solution
     CpGArray: An array containing all the possible CpG locations when maximizing the CpG locations
     If the entire array is considered, the most optimal solution may skip the required location.
    Instead, we will consider each subproblems where the first and last element are required locations
    and the elements in the middle are not required. Then the solutions can be concatenated.
     Outputs an array with the optimized CpG locations. This array will contain the required locations.
    */
    e.preventDefault();
    var count, end, i, lastNumber, locations, nec, necessaryCpGs, nums, opt, origLocation, prev, sol, solution, start, subArray;

    necessaryCpGs = JSON.parse(JSON.stringify(requiredCpG));

    nums = CpGArray.length;
    lastNumber = [CpGArray[nums - 1] + 13];
    locations = [0];
  
    for (var i = 0; i < CpGArray.length; i+=1){
        locations.push(CpGArray[i]);
    }
    locations.push(lastNumber);

    nums = locations.length;
    necessaryCpGs.reverse();
    start = 1;
    solution = [];
    count = 0;
  
    for (var k = 0; k < locations.length + 1; k += 1) {
        if (!!necessaryCpGs.length) {
            nec = necessaryCpGs.pop();
            end = locations.indexOf(nec) + 1;
        } else {
            end = locations.length;
        }

        subArray = locations.slice(start - 1, end);
        opt = [];
        prev = [];
        sol = [];
        opt.push(0);
        prev.push(1);
  
        for (var i = 1; i < subArray.length; i += 1) {
            opt.push(99999999999999999999999999999);
            prev.push(0);
  
            for (var j = 0; j < i+1; j += 1) {
                if (opt[i] > opt[j] + Math.pow(13 - locations[count + i] + locations[count + j], 2)) {
                    opt[i] = opt[j] + Math.pow(13 - locations[count + i] + locations[count + j], 2);
                    prev[i] = j;
                }
            }
        }
  
        i = subArray.length - 1;
        origLocation = locations.slice(start, end);
  
        while (i > 0) {
            sol.push(subArray[i]);
            i = prev[i];
        }
  
        sol.reverse();
  
        for (var i = 0; i < sol.length; i += 1){
            if (sol[i] > 0 && sol[i] < lastNumber) {
                solution.push(sol[i]);
            }
        }
  
        count += origLocation.length;
        start = end;
    }
    return solution;
}

function optimizeSequence(e, proteinSeq, solution) {
    /*
    This function will has two inputs:
          
          proteinSeq: String containing the original viruses protein sequence
          
          solution: An array with the optimal locations of the CpGs (output from optimizeSpacing)
      
    Outputs a string of DNA sequences that contains CpGs at the optimized CpG locations.
    */
    e.preventDefault();
    proteinIndex = Array.apply(null, Array(solution.length)).map(function (x, i) { return 0; })
    proteinZZ = Array.apply(null, Array(solution.length)).map(function (x, i) { return 0; })
    optimizedSequence = ""
    req = 0
  
    for (var i = 0; i < solution.length; i ++) {
        proteinIndex[i] = Math.floor(solution[i] / 3);
        proteinZZ[i] = solution[i] % 3;
    }
    x = proteinZZ.length - 1;
  
    proteinIndex.reverse();
    proteinZZ.reverse();
      
    followingAA = -1;
  
    for (var i = 0; i < proteinSeq.length; i ++) {
        followingAA += 1;
        a = proteinZZ[x];
      
        if (proteinIndex.includes(i) === false && proteinIndex.includes(followingAA+1) === false) {
            optimizedSequence += aminoToNuc_AT[proteinSeq[i]];
  
        } else if (proteinIndex.includes(i) === false && proteinIndex.includes(followingAA+1) === true) {
            if (a === 0){
                if (aminoToNuc_C[proteinSeq[i]] === "CGG") {
                    optimizeSequence += "AGG";
                } else {
                    optimizedSequence += aminoToNuc_C[proteinSeq[i]];
                }
            } else {
                optimizedSequence += aminoToNuc_AT[proteinSeq[i]];
            }
  
        } else if (proteinIndex.includes(i) === true && proteinIndex.includes(followingAA+1) === false) {
            if (a === 1 || a === 2){
                optimizedSequence += aminoToNuc_CG[proteinSeq[i]];
            } else {
                if (aminoToNuc_AT[proteinSeq[i]] === "AGG") {
                    optimizeSequence += "CGA";
                } else {
                    optimizedSequence += aminoToNuc_AT[proteinSeq[i]];
                }
            }
            x -= 1;

  
        } else if (proteinIndex.includes(i) === true && proteinIndex.includes(followingAA+1) === true) {
            if (proteinZZ[x-1] == 0){
                optimizedSequence += aminoToNuc_C[proteinSeq[i]];
            } else {
                optimizedSequence += aminoToNuc_CG[proteinSeq[i]];
            }
            x -= 1;
        }
    }
    return optimizedSequence;
}

function averageDistance(e, array) {
    /*
    Takes an array of locations as an input
    Outputs the average distance between each
    */
    e.preventDefault();
    var average, distance, sum;

    if (array.length == 0) {
        return;
    }

    var distance = Array.apply(null, Array(array.length-1)).map(function (x, i) { return 0; })
    sum = 0;
  
    for (var i = 0, lenArray = array.length - 1; i < lenArray; i += 1) {
        distance[i] = Math.abs(array[i + 1] - array[i]);
        sum += distance[i];
    }
  
    average = sum / distance.length;
    return average;
}