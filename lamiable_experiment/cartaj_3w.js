/* converts (uppercase) nucleotide letters to integers */
function n2i(n) {
    switch(n) {
    case 'A':
	return 2;
    case 'C':
	return 3;
    case 'G':
	return 4;
    case 'U':
	return 5;
    default:
	return -1;
    }
}

/* converts integers to (uppercase) nucleotide letters*/
function i2n(i) {
    switch(i) {
    case 2:
	return 'A';
    case 3:
	return 'C';
    case 4:
	return 'G';
    case 5:
	return 'U';
    default:
	return -1;
    }
}

function Strand(strand_start)
{
    this.strand_start = strand_start;
    this.before = [];
    this.unpaired = [];
    this.after = [];
    this.ids = [];

    this.get = function(i) {
	if (i < this.before.length)
	    return this.before[i];
	else if (i < this.before.length + this.unpaired.length)
	    return this.unpaired[i - this.before.length];
	else return this.after[i - this.before.length - this.unpaired.length];
    }

    this.size = function() {
	return this.before.length + this.unpaired.length + this.after.length;
    }
}

function Junction(size)
{
    this.size = size;
    this.strands = [];
    for(var i=0 ; i<size ; i++) {
	this.strands[i] = new Strand(-1);
    }

    this.last_WC5 = function (i) {
	var n = this.strands[i].before.length;
	return this.strands[i].before[n-1];
    };

    this.last_WC3 = function (i) {
	return this.strands[i].after[0];
    };

    this.n_nucleotides = function() {
	s = 0;
	for(i=0 ; i<this.size ; i++) s += this.strands[i].size();
	return s;
    }

    this.opposite_strand = function(i) {
	if (i == 0) {
	    return this.strands[this.strands.length-1];
	} else {
	    return this.strands[i-1];
	}
    }
}

function junction_from_sequence(seq)
{
    var a = seq.split("\n");
    var j = new Junction(a.length);
    var str;
    for(var i=0 ; i<a.length ; i++) {
	str = a[i].split("\.");
	j.strands[i].before = str[0].split('');
	j.strands[i].unpaired = str[1].split('');
	j.strands[i].after = str[2].split('');
    }
    return j;
}

function StackingTest()
{
    this.name = 'Stacking';
    this.compute = function (j, f, r) {
	var a = j.strands[0].unpaired.length;
	var b = j.strands[2].unpaired.length;
	var c = j.strands[1].unpaired.length;

        var sa = 1.0/(1+a);
        var sb = 1.0/(1+b);
        var sc = 1.0/(1+c);

	(a <= b) ? sb += 1.0/(1+b-a) : sa += 1.0/(1+a-b);
	(b <= c) ? sc += 1.0/(1+c-b) : sb += 1.0/(1+b-c);
	(c <= a) ? sa += 1.0/(1+a-c) : sc += 1.0/(1+c-a);

	(a <= b) ? sa++ : sb++;
	(b <= c) ? sb++ : sc++;
	(c <= a) ? sc++ : sa++;

	if (a == 0) sa = 10.0;
	if (b == 0) sb = 10.0;
	if (c == 0) sc = 10.0;

	switch(r) {
	case 0:
	    return sa;
	case 1:
	    return sb;
	case 2:
	    return sc;
	}
    }
}

function LengthTest()
{
    this.name = 'Lengths';
    this.compute = function (j, f, r) {
	var delta = (j.strands[(4-r)%3].unpaired.length -
		     j.strands[(5-r)%3].unpaired.length);

	switch(f) {
	case 'A':
	    break;
	case 'B':
	    delta = Math.abs(delta);
	    break;
	case 'C':
	    delta = -delta;
	    break;
	}

	switch(f) {
	case 'A':
	    break;
	case 'B':
	    switch(delta) {
	    case 0:
	    case 1:
		delta = 5.0;
		break;
	    case 2:
		delta = 1.0;
		break;
	    default:
		delta = 3.0 - delta;
		break;
	    }
	    break;
	case 'C':
	    extra = j.strands[(5-r)%3].unpaired.length - 5
	    if (extra > 0) delta += extra;
	    break;
	}

	if (delta > 10.0)
	    return 10.0;
	else if (delta < -10.0)
	    return -10.0;
	else return delta;
    }
}

function BonusTest()
{
    this.name = 'Bonus';
    this.compute = function (j, f, r) {
        bonus = 0.0;
        iV = (3-r)%3;
        iR = (4-r)%3;
        iB = (5-r)%3;

	switch(f) {
	case 'A':
	    // bonus 1
	    last_wc5 = j.strands[iR].before.length - 1;
	    start = last_wc5 + 2;
	    end = last_wc5 + 5;
	    for(var i=start ; i<=end ; i++) {
		if (j.strands[iR].get(i) == 'A') {
		    bonus += 1;
		    break;
		}
	    }

	    // bonus 2
	    first_non_wc_v3 = j.strands[iV].before.length + j.strands[iV].unpaired.length -1;
	    first_non_wc_r5 = j.strands[iR].before.length -1;
	    if (j.strands[iV].get(first_non_wc_v3) == 'A') {
		bonus++;
		if (j.strands[iR][first_non_wc_r5] == 'G') {
		    bonus += 0.5;
		}
	    }
	    break;
	case 'B':
	    bonus = 0.0;
	    break;
	case 'C':
	    first_non_wc_b3 = j.strands[iB].before.length + j.strands[iB].unpaired.length -1;
	    a = u = 0;
	    for(var i=first_non_wc_b3-1 ; i>=first_non_wc_b3-4 && i>=0 ; i--) {
		switch(j.strands[iB].get(i)) {
		case 'A':
		    a++;
		    break;
		case 'U':
		    u++;
		    break;
		}
	    }
	    if (a > 0)
		bonus = (a+u)/2.0;
	    else
		bonus = 0.0;
	    break;
	}

	return bonus;
    }
}

function ClosingTest()
{
    this.get_row = function(a, b) {
	if (a+b == 7)
	    return a-2;
	else if (i2n(a) == 'G' && i2n(b) == 'U')
	    return 4;
	else if (i2n(a) == 'U' && i2n(b) == 'G')
	    return 5;
	else throw 'error';
    }

    this.name = 'Closing BPs';
    this.compute = function (j, fam, r) {
	b5r3 = [[-0.08, -0.2, -0.2],
		[0.52, 0.485714285714286, 0.04],
		[0.16, 0.142857142857143, 0.2],
		[-0.2, -0.0285714285714286, 0.28],
		[-0.2, -0.2, -0.12], [-0.2, -0.2, -0.2]];

	r5v3 = [[-0.2, -0.2, -0.04],
		[0.4, 0.485714285714286, -0.04],
		[0.28, 0.142857142857143, 0.36],
		[-0.08, -0.0285714285714286, -0.12],
		[-0.2, -0.2, -0.2],
		[-0.2, -0.2, 0.04]];

	v5b3 = [[0.04, -0.0285714285714286, -0.04],
		[0.64, 0.314285714285714, 0.44],
		[-0.2, 0.142857142857143, -0.04],
		[-0.2, -0.0285714285714286, 0.04],
		[-0.2, -0.2, -0.2],
		[-0.08, -0.2, -0.2]]

	var f;
	switch(fam) {
	case 'A':
	    f = 0;
	    break;
	case 'B':
	    f = 1;
	    break;
	case 'C':
	    f = 2;
	    break;
	default:
	    throw 'error';
	    break;
	}

        iV = (3-r)%3;
        iR = (4-r)%3;
        iB = (5-r)%3;

	v5 = n2i(junction.last_WC5(iV));
        r5 = n2i(junction.last_WC5(iR));
        b5 = n2i(junction.last_WC5(iB));
	
        v3 = n2i(junction.last_WC3(iV));
        r3 = n2i(junction.last_WC3(iR));
        b3 = n2i(junction.last_WC3(iB));

	return (r5v3[this.get_row(r5,v3)][f] +
		b5r3[this.get_row(b5,r3)][f] +
		v5b3[this.get_row(v5,b3)][f]);
    }
}

function Eval()
{
    this.tests = [];
    this.coefs = [];
    this.values = [];
    this.totals = [];

    this.families = ['A', 'B', 'C'];
    this.rotations = [0, 1, 2];

    this.add_test = function(t, c) {
	this.tests.push(t);
	this.coefs.push(c);
    }

    this.compute = function(j) {
	this.values.length = 0;
	this.totals = [ [0, 0, 0], [0, 0, 0], [0, 0, 0] ];

	// computes the partial scores for all tests
	for(var i=0 ; i<this.tests.length ; i++) {
	    this.values[i] = [ [[], [], []],
			       [[], [], []],
			       [[], [], []] ];
	    for(var f=0 ; f<this.families.length ; f++)
		for(var r=0 ; r<this.rotations.length ; r++)
		    this.values[i][f][r] = this.tests[i].compute(j,this.families[f],r);
	}

	// computes the weighted sum and divides by the sum of the ranks
	var ranks;
	for(var f=0 ; f<this.families.length ; f++) {
	    for(var r=0 ; r<this.rotations.length ; r++) {
		ranks = 0;
		for(var i=0 ; i<this.tests.length ; i++) {
		    this.totals[f][r] += this.values[i][f][r] * this.coefs[i];
		    ranks += this.rank(i,f,r);
		}
		this.totals[f][r] /= ranks;
	    }
	}
    }

    // rank of configuration (fam,rot) for test i
    this.rank = function(i, fam, rot) {
	var rank = 1;
	var val = this.values[i][fam][rot];
	for(var f=0 ; f<this.families.length ; f++) {
	    for(var r=0 ; r<this.rotations.length ; r++) {
		if (this.values[i][f][r] > val) {
		    rank++;
		}
	    }
	}
	return rank;
    }

    // outputs an HTML table with all the scoring data
    this.toHtml = function() {
	var str = "<table>\n<tr>\n<th>Family</th><th>Rotation</th>";
	for(var i=0 ; i<this.tests.length ; i++) {
	    str += "<th>"+this.tests[i].name + "</th>";
	}
	str += "<th>Total</th></tr>\n";
	for(var r=0 ; r<this.rotations.length ; r++) {
	    for(var f=0 ; f<this.families.length ; f++) {
		str += "<tr onMouseOver=\"mydraw("+r+",'"+this.families[f]+"')\">\n";
		str += "<td>" + this.families[f] + "</td>\n";
		str += "<td>" + this.rotations[r] + "</td>\n";
		for(var i=0 ; i<this.tests.length ; i++) {
		    str += "<td>" + (new Number(this.values[i][f][r])).toFixed(2) + "</td>\n";
		}
		str += "<td>"+(new Number(this.totals[f][r])).toFixed(5) + "</td>";
		str += "</tr>\n";
	    }
	}
	str += "</table>";
	str += "<p>Move your mouse over the lines to see the junction drawn in the corresponding configuration.</p>";
	return str;
    }
}