(function (exports) {
    'use strict';

    exports.VERSION = "1.0.0-SNAPSHOT";

    var MU = 398600.5;

    var SECONDS_PER_DAY = 86400;

    var TWO_PI = 2 * Math.PI;

    var EQUATORIAL_RADIUS = 6378.137;

    var POLAR_RADIUS = 6356.7523142;

    var RAD2DEG = 180 / Math.PI;

    var DEG2RAD = Math.PI / 180;

    Object.defineProperty(Function.prototype, '__doc__', {
        get: function () {
            var comment = this.toString(),
                __doc__ = '';
            if (comment = comment.match(/\/\*[!*]([\s\S]*?)\*\//)) {
                __doc__ = comment[1];
                __doc__ = __doc__.replace(/\*/g, ' ').replace(/\ +/g, ' ');
            }
            return __doc__.trim();
        }
    });

    var copyArray = function (a) {
        var o = new Array(a.length);
        for (var i = 0; i < a.length; i++) {
            o[i] = a[i];
        }
        return o;
    }

    var magnitude = function (a) {
        return Math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    }

    var normalize = function (a) {
        var m = magnitude(a);
        return [a[0] / m, a[1] / m, a[2] / m];
    }

    var scale = function (a, s) {
        return [a[0] * s, a[1] * s, a[2] * s];
    }

    var gravity = function (a) {
        var v = normalize(a);
        var m = magnitude(a);
        var g = -(MU / (m * m));
        return scale(v, g);
    }

    var julainDate = function (t) {
        return (t / (SECONDS_PER_DAY * 1000)) + 2440587.5;
    }

    var j2kDays = function (t) {
        return julainDate(t) - 2451545;
    }

    var gmst = function (t) {
        var D = j2kDays(t);
        var theta = (4.894961213 + 6.300388099 * D) % TWO_PI;
        return (theta < 0) ? TWO_PI + theta : theta;
    }

    exports.now = function () {
        /**
         * () -> number
         * 
         * Return the current unix timestamp (milliseconds).
         */
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
        /**
         * (number, [number, number, number]) -> [number, number, number]
         * 
         * Convert ECI coordinates (kilometers) to ECEF coordinates
         * (kilometers), at the given unix time (milliseconds).
         */
        var a = gmst(t);
        return rotateZ(-a, eci);
    }

    exports.ecef2eci = function (t, ecef) {
        /**
         * (number, [number, number, number]) -> [number, number, number]
         * 
         * Convert ECEF coordinates (kilometers) to ECI coordinates
         * (kilometers), at the given unix time (milliseconds).
         */
        var a = gmst(t);
        return rotateZ(a, ecef);
    }

    exports.ecef2geodetic = function (ecef) {
        /**
         * ([number, number, number]) -> [number, number, number]
         * 
         * Convert ECEF coordinates (kilometers) to geodetic latitude (degrees),
         * longitude (degrees), and altitude (kilometers).
         */
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
        /**
         * ([number, number, number]) -> [number, number, number]
         * 
         * Convert geodetic latitude (degrees), longitude (degrees), and
         * altitude (kilometers) to ECEF coordinates (kilometers).
         */
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
        this._initEpoch = epoch;
        this._initPosition = position;
        this._initVelocity = velocity;
        this._step = step || 1;
        this.reset();
    }

    exports.Propagator.prototype.reset = function () {
        this._epoch = this._initEpoch;
        this._position = copyArray(this._initPosition);
        this._velocity = copyArray(this._initVelocity);
        return this;
    }

    exports.Propagator.prototype.setStep = function (seconds) {
        if (seconds > 0) {
            this._step = seconds;
        } else {
            throw new Error("Step size must be a positive number.");
        }
        return this;
    }

    exports.Propagator.prototype.propagate = function (epoch) {
        var dir = (epoch < this._epoch) ? -1 : 1;
        while (epoch !== this._epoch) {
            var s = this._position;
            var v = this._velocity;
            var a = gravity(s);
            var dist = Math.abs(epoch - this._epoch) / 1000;
            var t = Math.min(this._step, dist) * dir;
            this._epoch += t * 1000;
            var sNext = [
                s[0] + (v[0] * t) + (0.5 * a[0] * t * t),
                s[1] + (v[1] * t) + (0.5 * a[1] * t * t),
                s[2] + (v[2] * t) + (0.5 * a[2] * t * t)
            ];
            var vNext = [
                v[0] + (a[0] * t),
                v[1] + (a[1] * t),
                v[2] + (a[2] * t)
            ];
            this._position = sNext;
            this._velocity = vNext;
        }
        return this;
    }

    exports.Propagator.prototype.getEpoch = function () {
        return this._epoch;
    }

    exports.Propagator.prototype.getECI = function () {
        return this._position;
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
})(this.GeoHack = this.GeoHack || {});
