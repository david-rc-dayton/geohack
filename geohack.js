var GeoHack = {};

(function (exports) {
    'use strict';

    exports.VERSION = "1.0.0";

    var MU = 398600.5;

    var SECONDS_PER_DAY = 86400;

    var TWO_PI = 2 * Math.PI;

    var EQUATORIAL_RADIUS = 6378.137;

    var FLATTENING = 1 / 298.257223563;

    var POLAR_RADIUS = EQUATORIAL_RADIUS * (1 - FLATTENING);

    var SOLAR_DAY = 86400

    var SIDEREAL_DAY = 0.99726958 * SOLAR_DAY

    var RAD2DEG = 180 / Math.PI;

    var DEG2RAD = Math.PI / 180;

    var clone = function (a) {
        return a.slice(0);
    }

    var magnitude = function (a) {
        return Math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    }

    var normalize = function (a) {
        var m = magnitude(a);
        return [a[0] / m, a[1] / m, a[2] / m];
    }

    var scale = function (s, a) {
        return [a[0] * s, a[1] * s, a[2] * s];
    }

    var add = function () {
        var state = clone(arguments[0]);
        for (var i = 1; i < arguments.length; i++) {
            for (var j = 0; j < state.length; j++) {
                state[j] += arguments[i][j];
            }
        }
        return state;
    }

    var gravity = function (a) {
        var n = normalize(a);
        var r = magnitude(a);
        var g = -(MU / (r * r));
        return scale(g, n);
    }

    var julianDate = function (t) {
        return (t / (SECONDS_PER_DAY * 1000)) + 2440587.5;
    }

    var j2kDays = function (t) {
        return julianDate(t) - 2451545;
    }

    var gmst = function (t) {
        var D = j2kDays(t);
        var theta = (4.894961213 + 6.300388099 * D) % TWO_PI;
        return (theta < 0) ? TWO_PI + theta : theta;
    }

    exports.now = function () {
        return new Date().getTime();
    }

    var rotateZ = function (theta, v) {
        return [
            v[0] * Math.cos(theta) + v[1] * -Math.sin(theta) + v[2] * 0,
            v[0] * Math.sin(theta) + v[1] * Math.cos(theta) + v[2] * 0,
            v[0] * 0 + v[1] * 0 + v[2] * 1
        ];
    }

    exports.eci2ecef = function (t, eci) {
        var a = gmst(t);
        return rotateZ(-a, eci);
    }

    exports.ecef2eci = function (t, ecef) {
        var a = gmst(t);
        return rotateZ(a, ecef);
    }

    exports.ecef2geodetic = function (ecef) {
        var x = ecef[0];
        var y = ecef[1];
        var z = ecef[2];
        var a = EQUATORIAL_RADIUS;
        var b = POLAR_RADIUS;
        var R = Math.sqrt((x * x) + (y * y));
        var f = (a - b) / a;
        var e2 = (2 * f) - (f * f);
        var lon = Math.atan2(y, x);
        var lat = Math.atan2(z, R);
        var c = 0;
        for (var k = 0; k < 20; k++) {
            c = 1 / Math.sqrt(1 - (e2 * (Math.sin(lat) * Math.sin(lat))));
            lat = Math.atan2(z + (a * c * e2 * Math.sin(lat)), R);
        }
        var alt = (R / Math.cos(lat)) - (a * c);
        return [lat * RAD2DEG, lon * RAD2DEG, alt];
    }

    exports.geodetic2ecef = function (geodetic) {
        var lat = geodetic[0] * DEG2RAD;
        var lon = geodetic[1] * DEG2RAD;
        var alt = geodetic[2];
        var a = EQUATORIAL_RADIUS;
        var b = POLAR_RADIUS;
        var f = (a - b) / a;
        var e2 = (2 * f) - (f * f);
        var n = a / Math.sqrt(1 - (e2 * (Math.sin(lat) * Math.sin(lat))));
        var x = (n + alt) * Math.cos(lat) * Math.cos(lon);
        var y = (n + alt) * Math.cos(lat) * Math.sin(lon);
        var z = ((n * (1 - e2)) + alt) * Math.sin(lat);
        return [x, y, z];
    }

    exports.topocentric = function (geoOrigin, ecefTarget) {
        var lat = geoOrigin[0] * DEG2RAD;
        var lon = geoOrigin[1] * DEG2RAD;
        var oECEF = exports.geodetic2ecef(geoOrigin);
        var tECEF = ecefTarget;
        var rX = tECEF[0] - oECEF[0];
        var rY = tECEF[1] - oECEF[1];
        var rZ = tECEF[2] - oECEF[2];
        var S = (
            (Math.sin(lat) * Math.cos(lon) * rX) +
            (Math.sin(lat) * Math.sin(lon) * rY) -
            (Math.cos(lat) * rZ)
        );
        var E = (-Math.sin(lon) * rX) + (Math.cos(lon) * rY);
        var Z = (
            (Math.cos(lat) * Math.cos(lon) * rX) +
            (Math.cos(lat) * Math.sin(lon) * rY) +
            (Math.sin(lat) * rZ)
        );
        return [S, E, Z];
    }

    exports.lookangle = function (geoOrigin, ecefTarget) {
        var sez = exports.topocentric(geoOrigin, ecefTarget);
        var s = sez[0];
        var e = sez[1];
        var z = sez[2];
        var dist = Math.sqrt((s * s) + (e * e) + (z * z));
        var el = Math.asin(z / dist);
        var az = Math.atan2(-e, s) + Math.PI;
        return [az * RAD2DEG, el * RAD2DEG, dist];
    }

    exports.Propagator = function (epoch, position, velocity, step) {
        this.initEpoch = epoch;
        this.initPosition = position;
        this.initVelocity = velocity;
        this.step = step || 10;
        this.reset();
    }

    exports.Propagator.prototype.reset = function () {
        this.epoch = this.initEpoch;
        this.position = clone(this.initPosition);
        this.velocity = clone(this.initVelocity);
        return this;
    }

    exports.Propagator.prototype.setStep = function (seconds) {
        if (seconds > 0) {
            this.step = seconds;
        } else {
            throw new Error("Step size must be a positive number.");
        }
        return this;
    }

    exports.Propagator.prototype.propagate = function (newEpoch) {
        while (this.epoch !== newEpoch) {
            var direction = (newEpoch > this.epoch) ? 1 : -1;
            var delta = Math.abs(newEpoch - this.epoch) / 1000;
            var r = this.position;
            var v = this.velocity;
            var h = Math.min(this.step, delta) * direction;
            var kv1 = gravity(r);
            var kr1 = v;
            var kv2 = gravity(add(r, scale(h / 2, kr1)));
            var kr2 = add(v, scale(h / 2, kv1));
            var kv3 = gravity(add(r, scale(h / 2, kr2)));
            var kr3 = add(v, scale(h / 2, kv2));
            var kv4 = gravity(add(r, scale(h, kr3)));
            var kr4 = add(v, scale(h, kv3));
            var vpost = add(kv1, scale(2, kv2), scale(2, kv3), kv4);
            var rpost = add(kr1, scale(2, kr2), scale(2, kr3), kr4);
            this.velocity = add(v, scale(h / 6, vpost));
            this.position = add(r, scale(h / 6, rpost));
            this.epoch += h * 1000;
        }
        return this;
    }

    exports.Propagator.prototype.getEpoch = function () {
        return this.epoch;
    }

    exports.Propagator.prototype.getECI = function () {
        return this.position;
    }

    exports.Propagator.prototype.getECEF = function () {
        return exports.eci2ecef(this.getEpoch(), this.getECI());
    }

    exports.Propagator.prototype.getGeodetic = function () {
        return exports.ecef2geodetic(this.getECEF());
    }

    exports.Propagator.prototype.getGeodeticDegrees = function () {
        var c = this.getGeodetic();
        return [c[0] * RAD2DEG, c[1] * RAD2DEG, c[2]];
    }

    exports.Propagator.prototype.getLookAngle = function (geoOriginDeg) {
        var ecef = this.getECEF();
        return exports.lookangle(geoOriginDeg, ecef);
    }

    exports.testSatellite = new exports.Propagator(
        exports.now(),
        [42164.17207, 0, 0],
        [0, 3.074660235, 0]
    );
})(GeoHack);

var issEpoch = Date.UTC(2017, 0, 14, 20, 55);
var issPosition = [-111.31690, 6662.17019, 1255.34162];
var issVelocity = [-4.813744845, -1.196254215, 5.853086769];
var issPropagator = new GeoHack.Propagator(issEpoch, issPosition, issVelocity);
var issNextEpoch = Date.UTC(2017, 0, 15, 20, 55);
var issNextPosition = [945.92445, -6089.61961, -2849.12060];
var issTestPosition = issPropagator.propagate(issNextEpoch).position;
console.log(issNextPosition);
console.log(issTestPosition);