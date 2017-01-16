var GeoHack = require('./geohack');
var assert = require('assert');

var LIMIT = 5 // kilometers

function distance(a, b) {
    return Math.sqrt(
        (a[0] - b[0]) * (a[0] - b[0]) +
        (a[1] - b[1]) * (a[1] - b[1]) +
        (a[2] - b[2]) * (a[2] - b[2])
    );
}

function fix(n) {
    return n.toFixed(3);
}

var issEpoch = Date.UTC(2017, 0, 14, 20, 55);
var issPosition = [-111.31690, 6662.17019, 1255.34162];
var issVelocity = [-4.813744845, -1.196254215, 5.853086769];
var issPropagator = new GeoHack.Propagator(issEpoch, issPosition, issVelocity);
var issNextEpoch = Date.UTC(2017, 0, 15, 20, 55);
var issNextPosition = [945.92445, -6089.61961, -2849.12060];
var issTestPosition = issPropagator.propagate(issNextEpoch).position;
var deltaPosition = distance(issNextPosition, issTestPosition);
console.log("Expected:  ", issNextPosition.map(fix), "km");
console.log("Actual:    ", issTestPosition.map(fix), "km");
console.log("Difference:", fix(deltaPosition), "km");
assert(deltaPosition <= LIMIT);
