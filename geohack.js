var GeoHack = GeoHack || {};

// Constants //////////////////////////////////////////////////////////////////

GeoHack.VERSION = "1.0.0-SNAPSHOT";

GeoHack.MU = 398600.5;

GeoHack.SECONDS_PER_DAY = 86400;

GeoHack.TWO_PI = 2 * Math.PI;

GeoHack.EQUATORIAL_RADIUS = 6378.137;

GeoHack.POLAR_RADIUS = 6356.7523142;

GeoHack.RAD2DEG = 180 / Math.PI;

GeoHack.DEG2RAD = Math.PI / 180;

// Utilities //////////////////////////////////////////////////////////////////

GeoHack.copyArray = function (a) {
    var o = new Array(a.length);
    for (var i = 0; i < a.length; i++) {
        o[i] = a[i];
    }
    return o;
}

GeoHack.magnitude = function (a) {
    return Math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

GeoHack.normalize = function (a) {
    var m = GeoHack.magnitude(a);
    return [a[0] / m, a[1] / m, a[2] / m];
}

GeoHack.scale = function (a, s) {
    return [a[0] * s, a[1] * s, a[2] * s];
}

GeoHack.gravity = function (a) {
    var v = GeoHack.normalize(a);
    var m = GeoHack.magnitude(a);
    var g = -(GeoHack.MU / (m * m));
    return GeoHack.scale(v, g);
}

GeoHack.julainDate = function (t) {
    return (t / (GeoHack.SECONDS_PER_DAY * 1000)) + 2440587.5;
}

GeoHack.j2kDays = function (t) {
    return GeoHack.julainDate(t) - 2451545;
}

GeoHack.gmst = function (t) {
    var D = GeoHack.j2kDays(t);
    var theta = (4.894961213 + 6.300388099 * D) % GeoHack.TWO_PI;
    return (theta < 0) ? GeoHack.TWO_PI + theta : theta;
}

GeoHack.rotateZ = function (theta, v) {
    return [
        v[0] * Math.cos(theta) + v[1] * -Math.sin(theta) + v[2] * 0,
        v[0] * Math.sin(theta) + v[1] * Math.cos(theta) + v[2] * 0,
        v[0] * 0 + v[1] * 0 + v[2] * 1
    ];
}

// Coordinates ////////////////////////////////////////////////////////////////

GeoHack.eci2ecef = function (t, v) {
    var a = GeoHack.gmst(t);
    return GeoHack.rotateZ(-a, v);
}

GeoHack.ecef2geodetic = function (v) {
    var x = v[0];
    var y = v[1];
    var z = v[2];
    var a = GeoHack.EQUATORIAL_RADIUS;
    var b = GeoHack.POLAR_RADIUS;
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
    return [lat, lon, alt];
}

GeoHack.geodetic2ecef = function (v) {
    var lat = v[0];
    var lon = v[1];
    var alt = v[2];
    var a = GeoHack.EQUATORIAL_RADIUS;
    var b = GeoHack.POLAR_RADIUS;
    var f = (a - b) / a;
    var e2 = (2 * f) - (f * f);
    var n = a / Math.sqrt(1 - (e2 * (Math.sin(lat) * Math.sin(lat))));
    var x = (n + alt) * Math.cos(lat) * Math.cos(lon);
    var y = (n + alt) * Math.cos(lat) * Math.sin(lon);
    var z = ((n * (1 - e2)) + alt) * Math.sin(lat);
    return [x, y, z];
}

GeoHack.topocentric = function (geoOrigin, ecefTarget) {
    var lat = geoOrigin[0];
    var lon = geoOrigin[1];
    var oECEF = GeoHack.geodetic2ecef(geoOrigin);
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

GeoHack.lookAngle = function (geoOrigin, ecefTarget) {
    var sez = GeoHack.topocentric(geoOrigin, ecefTarget);
    var s = sez[0];
    var e = sez[1];
    var z = sez[2];
    var dist = Math.sqrt((s * s) + (e * e) + (z * z));
    var el = Math.asin(z / dist);
    var az = Math.atan2(-e, s) + Math.PI;
    return [az, el, dist];
}

// Propagator /////////////////////////////////////////////////////////////////

GeoHack.Propagator = function (epoch, position, velocity, step) {
    this._initEpoch = epoch;
    this._initPosition = position;
    this._initVelocity = velocity;
    this._step = step || 0.1;
    this.reset();
}

GeoHack.Propagator.prototype.reset = function () {
    this._epoch = this._initEpoch;
    this._position = GeoHack.copyArray(this._initPosition);
    this._velocity = GeoHack.copyArray(this._initVelocity);
    return this;
}

GeoHack.Propagator.prototype.setStep = function (seconds) {
    if (seconds > 0) {
        this._step = seconds;
    } else {
        throw new Error("Step size must be a positive number.");
    }
    return this;
}

GeoHack.Propagator.prototype.propagate = function (epoch) {
    var dir = (epoch < this._epoch) ? -1 : 1;
    while (epoch !== this._epoch) {
        var s = this._position;
        var v = this._velocity;
        var a = GeoHack.gravity(s);
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

GeoHack.Propagator.prototype.getEpoch = function () {
    return this._epoch;
}

GeoHack.Propagator.prototype.getECI = function () {
    return this._position;
}

GeoHack.Propagator.prototype.getECEF = function () {
    return GeoHack.eci2ecef(this.getEpoch(), this.getECI());
}

GeoHack.Propagator.prototype.getGeodetic = function () {
    return GeoHack.ecef2geodetic(this.getECEF());
}

GeoHack.Propagator.prototype.getGeodeticDegrees = function () {
    var c = this.getGeodetic();
    return [c[0] * GeoHack.RAD2DEG, c[1] * GeoHack.RAD2DEG, c[2]];
}

// Test ///////////////////////////////////////////////////////////////////////

GeoHack.testEpochStart = new Date().getTime();

GeoHack.testEpochEnd = GeoHack.testEpochStart + 86164.09171 * 1000;

GeoHack.testProp = new GeoHack.Propagator(
    GeoHack.testEpochStart, [42164.17207, 0, 0], [0, 3.074660235, 0]);
