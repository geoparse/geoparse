(function () {
    var mapObj = {{ this.uid }};
    var defaultRadiusPx = {{ this.default_radius_px }};
    var radiusM = 0;
    var maxRadiusM = {{ this.max_radius_km }} * 1000;
    var minRadiusM = 50;
    // { alias: [fieldName, func] }
    //   func in {"sum","mean"/"avg","max","min","count","unique"}
    var aggDict = {{ this.agg_dict_json }};
    var aggAliases = Object.keys(aggDict);
    var circleOverlay, centerDot, badge, radiusLabel, initialized = false;
    var isResizing = false;
    var isMoving = false;
    var moveOffset = null;

    function haversineMeters(lat1, lon1, lat2, lon2) {
        var R = 6371000;
        function toRad(d) { return d * Math.PI / 180; }
        var dLat = toRad(lat2 - lat1), dLon = toRad(lon2 - lon1);
        var a = Math.sin(dLat / 2) ** 2 +
                Math.cos(toRad(lat1)) * Math.cos(toRad(lat2)) * Math.sin(dLon / 2) ** 2;
        return 2 * R * Math.asin(Math.sqrt(a));
    }

    // Folium/Leaflet's addData() always sets `.feature` on every marker it
    // creates from a GeoJson layer, so this is the one place that matters.
    function getProps(layer) {
        return (layer && layer.feature && layer.feature.properties) || {};
    }

    // Coerce a property value to a number, stripping thousands separators
    // if present. Returns null for anything that shouldn't contribute to a
    // numeric aggregation (blank, boolean, non-numeric string).
    function toNumber(v) {
        if (v === null || v === undefined || v === "") return null;
        if (typeof v === "number") return isNaN(v) ? null : v;
        if (typeof v === "boolean") return null;
        if (typeof v === "string") {
            var cleaned = v.replace(/,/g, "").trim();
            if (cleaned === "") return null;
            var n = Number(cleaned);
            return isNaN(n) ? null : n;
        }
        return null;
    }

    function formatNumber(n) {
        if (n === null || n === undefined || !isFinite(n)) return "-";
        if (Math.abs(n) >= 1000) {
            return n.toLocaleString(undefined, { maximumFractionDigits: 2 });
        }
        return String(Math.round(n * 100) / 100);
    }

    // Convert a screen-pixel distance at `center` into metres, at the
    // current zoom / latitude of the map.
    function pixelsToMeters(center, px) {
        var cp = mapObj.latLngToContainerPoint(center);
        var ep = L.point(cp.x + px, cp.y);
        var el = mapObj.containerPointToLatLng(ep);
        return haversineMeters(center.lat, center.lng, el.lat, el.lng);
    }

    // Work out where the circle should start: tucked into the top-left
    // corner, just below the measurement control. We measure the control's
    // actual bounding box so this works regardless of Leaflet version /
    // control size.
    function getStartCenter() {
        var container = mapObj.getContainer();
        var containerRect = container.getBoundingClientRect();

        // Extra vertical gap so the badge -- which hovers ~24 px above the
        // circle's top edge -- clears the measurement control instead of
        // covering it.
        var badgeClearance = 40;

        var measureEl = container.querySelector(".leaflet-control-measure")
                     || container.querySelector(".leaflet-control.leaflet-bar");
        var belowY, leftX;
        if (measureEl) {
            var r = measureEl.getBoundingClientRect();
            leftX = (r.left - containerRect.left) + r.width / 2;
            belowY = (r.bottom - containerRect.top) + defaultRadiusPx + badgeClearance;
        } else {
            // Fallback if the control isn't found: sit just inside the
            // top-left with the usual 10 px margin.
            leftX = 10 + defaultRadiusPx;
            belowY = 10 + 34 + defaultRadiusPx + badgeClearance;
        }

        leftX = Math.max(leftX, defaultRadiusPx + 4);
        belowY = Math.max(belowY, defaultRadiusPx + 4);

        return mapObj.containerPointToLatLng(L.point(leftX, belowY));
    }

    // Snap a raw metre value to the step that matches its magnitude. The
    // bucket is chosen from the raw value, so the step changes smoothly as
    // you cross thresholds.
    function snapRadius(rawM) {
        var step;
        if (rawM < 100)         step = 10;
        else if (rawM < 1000)   step = 50;
        else if (rawM < 10000)  step = 500;
        else if (rawM < 100000) step = 5000;
        else                    step = 50000;
        return Math.round(rawM / step) * step;
    }

    function formatRadius(m) {
        if (m < 1000) return Math.round(m) + " m";
        var km = m / 1000;
        if (km === Math.round(km)) return km.toFixed(0) + " km";
        return km.toFixed(1) + " km";
    }

    // Recursively walk every layer on the map (through FeatureGroups /
    // GeoJson / MarkerCluster) collecting point-like markers together with
    // their GeoJSON props. Polygons/lines never match the CircleMarker/
    // Marker checks, so buffers, choropleths, etc. are left out. The tool's
    // own circle and centre dot are excluded too. `seen` dedupes by Leaflet
    // layer id so a marker reached through more than one parent path is
    // only counted once.
    function collectPoints(layer, out, seen) {
        var id = L.Util.stamp(layer);
        if (seen.has(id)) return;
        seen.add(id);

        if ((layer instanceof L.CircleMarker || layer instanceof L.Marker) &&
            layer !== circleOverlay && layer !== centerDot) {
            out.push({ latlng: layer.getLatLng(), props: getProps(layer) });
        }
        if (typeof layer.getAllChildMarkers === "function") {
            layer.getAllChildMarkers().forEach(function (m) {
                var mid = L.Util.stamp(m);
                if (!seen.has(mid)) {
                    seen.add(mid);
                    out.push({ latlng: m.getLatLng(), props: getProps(m) });
                }
            });
            return;
        }
        if (typeof layer.eachLayer === "function") {
            layer.eachLayer(function (child) { collectPoints(child, out, seen); });
        }
    }

    // Reduce one aggregation entry's collected values according to func.
    function aggregateOne(values, func) {
        if (func === "count") {
            var n = 0;
            for (var i = 0; i < values.length; i++) {
                if (values[i] !== null && values[i] !== undefined && values[i] !== "") n++;
            }
            return n;
        }
        if (func === "unique") {
            var seen = new Set();
            for (var j = 0; j < values.length; j++) {
                var v = values[j];
                if (v === null || v === undefined || v === "") continue;
                seen.add(String(v));
            }
            return seen.size;
        }

        // Numeric-only aggregations: strip non-numeric values.
        var nums = [];
        for (var k = 0; k < values.length; k++) {
            var num = toNumber(values[k]);
            if (num !== null) nums.push(num);
        }
        if (!nums.length) return null;
        if (func === "sum") {
            var s = 0;
            for (var a = 0; a < nums.length; a++) s += nums[a];
            return s;
        }
        if (func === "mean" || func === "avg") {
            var tot = 0;
            for (var b = 0; b < nums.length; b++) tot += nums[b];
            return tot / nums.length;
        }
        if (func === "max") {
            var mx = nums[0];
            for (var c = 1; c < nums.length; c++) if (nums[c] > mx) mx = nums[c];
            return mx;
        }
        if (func === "min") {
            var mn = nums[0];
            for (var d = 1; d < nums.length; d++) if (nums[d] < mn) mn = nums[d];
            return mn;
        }
        return null;
    }

    // Counts points inside the circle and aggregates every entry in
    // aggDict over those same points.
    function summarizeCircle(center, radius) {
        var pts = [];
        var seen = new Set();
        mapObj.eachLayer(function (l) { collectPoints(l, pts, seen); });

        var count = 0;
        var buckets = {};
        aggAliases.forEach(function (alias) { buckets[alias] = []; });

        for (var i = 0; i < pts.length; i++) {
            if (haversineMeters(center.lat, center.lng,
                                pts[i].latlng.lat, pts[i].latlng.lng) <= radius) {
                count++;
                aggAliases.forEach(function (alias) {
                    var field = aggDict[alias][0];
                    buckets[alias].push(pts[i].props[field]);
                });
            }
        }

        var results = {};
        aggAliases.forEach(function (alias) {
            results[alias] = aggregateOne(buckets[alias], aggDict[alias][1]);
        });
        return { count: count, results: results };
    }

    // Circle radius in screen pixels, derived from its bounds.
    function getPixelRadius() {
        if (!circleOverlay) return 0;
        var b = circleOverlay.getBounds();
        var ne = mapObj.latLngToContainerPoint(b.getNorthEast());
        var sw = mapObj.latLngToContainerPoint(b.getSouthWest());
        return (ne.x - sw.x) / 2;
    }

    // True if the given container point is close enough to the circle's
    // stroke to count as "grabbing the border".
    function isNearBorder(containerPoint) {
        var c = circleOverlay.getLatLng();
        var centerPt = mapObj.latLngToContainerPoint(c);
        var dist = containerPoint.distanceTo(centerPt);
        var pr = getPixelRadius();
        var tolerance = Math.min(10, Math.max(3, pr * 0.4));
        return Math.abs(dist - pr) < tolerance;
    }

    // Refresh the badge (count + aggregations + position) and keep the
    // centre dot glued to the circle's centre.
    function updateBadge() {
        if (!circleOverlay) return;
        var c = circleOverlay.getLatLng();

        if (centerDot) centerDot.setLatLng(c);

        if (!badge) return;
        var summary = summarizeCircle(c, radiusM);
        var pt = mapObj.latLngToContainerPoint(c);
        var pr = getPixelRadius();

        badge.style.left = pt.x + "px";
        badge.style.top  = (pt.y - pr - 8) + "px";

        var lines = [summary.count + (summary.count === 1 ? " point" : " points")];
        aggAliases.forEach(function (alias) {
            lines.push(alias + ": " + formatNumber(summary.results[alias]));
        });
        badge.innerHTML = lines.join("<br>");
    }

    function showRadiusLabel(containerPoint) {
        if (!radiusLabel) return;
        radiusLabel.textContent = formatRadius(radiusM);
        radiusLabel.style.left = (containerPoint.x + 14) + "px";
        radiusLabel.style.top = (containerPoint.y + 14) + "px";
        radiusLabel.style.display = "block";
    }

    function hideRadiusLabel() {
        if (radiusLabel) radiusLabel.style.display = "none";
    }

    function setCursor(cursor) {
        if (circleOverlay && circleOverlay._path) {
            circleOverlay._path.style.cursor = cursor;
        }
    }

    function endDrag() {
        if (!isResizing && !isMoving) return;
        isResizing = false;
        isMoving = false;
        moveOffset = null;
        hideRadiusLabel();
        mapObj.dragging.enable();
        setCursor("move");
    }

    function initCircleTool() {
        if (initialized) return;
        initialized = true;

        var center = getStartCenter();
        radiusM = pixelsToMeters(center, defaultRadiusPx);

        circleOverlay = L.circle(center, {
            radius: radiusM, color: "#e53e3e", weight: 2,
            fillColor: "#feb2b2", fillOpacity: 0.2, interactive: true
        }).addTo(mapObj);

        // Non-interactive so clicks fall through to the circle's fill
        // (which drives the drag-to-move logic).
        centerDot = L.circleMarker(center, {
            radius: 4, color: "#ffffff", weight: 1.5,
            fillColor: "#e53e3e", fillOpacity: 1.0, interactive: false
        }).addTo(mapObj);

        badge = L.DomUtil.create("div", "", mapObj.getContainer());
        badge.style.cssText =
            "position:absolute; z-index:1000; background:rgba(255,255,255,0.95);" +
            "padding:3px 8px; border-radius:4px; font-family:sans-serif;" +
            "font-size:12px; font-weight:bold; color:#c53030; line-height:1.4;" +
            "pointer-events:none; white-space:nowrap;" +
            "box-shadow:0 1px 4px rgba(0,0,0,0.25);" +
            "transform:translate(-50%,-100%); border:1px solid #e53e3e;";
        badge.textContent = "0 points";

        radiusLabel = L.DomUtil.create("div", "", mapObj.getContainer());
        radiusLabel.style.cssText =
            "position:absolute; z-index:1001; display:none;" +
            "background:rgba(229,62,62,0.95); color:white;" +
            "padding:2px 7px; border-radius:4px;" +
            "font-family:sans-serif; font-size:12px; font-weight:bold;" +
            "pointer-events:none; white-space:nowrap;" +
            "box-shadow:0 1px 4px rgba(0,0,0,0.3);";
        radiusLabel.textContent = formatRadius(radiusM);

        circleOverlay.on("mouseover", function () {
            if (!isResizing && !isMoving) setCursor("move");
        });
        circleOverlay.on("mouseout", function () {
            if (!isResizing && !isMoving) setCursor("");
        });
        circleOverlay.on("mousemove", function (e) {
            if (isResizing || isMoving) return;
            setCursor(isNearBorder(e.containerPoint) ? "ew-resize" : "move");
        });

        circleOverlay.on("mousedown", function (e) {
            L.DomEvent.stopPropagation(e);
            var c = circleOverlay.getLatLng();

            if (isNearBorder(e.containerPoint)) {
                isResizing = true;
                setCursor("ew-resize");
                showRadiusLabel(e.containerPoint);
            } else {
                isMoving = true;
                setCursor("move");
                var centerPt = mapObj.latLngToContainerPoint(c);
                moveOffset = L.point(
                    centerPt.x - e.containerPoint.x,
                    centerPt.y - e.containerPoint.y
                );
            }
            mapObj.dragging.disable();
        });

        mapObj.on("mousemove", function (e) {
            if (isResizing) {
                var c = circleOverlay.getLatLng();
                var rawR = haversineMeters(c.lat, c.lng, e.latlng.lat, e.latlng.lng);
                rawR = Math.max(minRadiusM, Math.min(rawR, maxRadiusM));
                radiusM = snapRadius(rawR);
                circleOverlay.setRadius(radiusM);
                updateBadge();
                showRadiusLabel(e.containerPoint);
            } else if (isMoving) {
                var newPt = L.point(
                    e.containerPoint.x + moveOffset.x,
                    e.containerPoint.y + moveOffset.y
                );
                circleOverlay.setLatLng(mapObj.containerPointToLatLng(newPt));
                updateBadge();
            }
        });

        mapObj.on("mouseup", endDrag);
        document.addEventListener("mouseup", endDrag); // safety net for mouseup off-map

        // We listen to `zoomend`, NOT `zoom`: during the zoom animation the
        // map's projection has already jumped to the final zoom, but the
        // SVG overlay pane is still being CSS-scaled mid-flight, so pixel
        // math would be out of sync with what's on screen until it finishes.
        mapObj.on("move zoomend viewreset", updateBadge);

        updateBadge();
    }

    // Defer placing the circle until the map's initial fit_bounds (called
    // later, after _base_map returns, once layers/data are added) has
    // settled -- otherwise we'd grab Leaflet's pre-fit default center.
    // Falls back to a timeout for maps with nothing to fit.
    mapObj.once("moveend", function () { initCircleTool(); });
    setTimeout(function () { initCircleTool(); }, 800);
})();
