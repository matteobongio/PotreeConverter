
#pragma once

#include <algorithm>
#include <cstddef>
#include <execution>

#include "Vector3.h"
#include "structures.h"
#include "Attributes.h"
#include "PotreeConverter.h"


struct SamplerCurvature : public Sampler {

	size_t maxNeighbours;
	size_t samples;

	SamplerCurvature(size_t maxNeighbours, size_t samples) 
	:
		maxNeighbours(maxNeighbours),
		samples(samples)
	{
	}

	struct HeapCmp {
		const vector<double>& dist;
		HeapCmp(const vector<double>& dist) : dist(dist) {}
		bool operator()(size_t a, size_t b) const {
			return dist[a] < dist[b];
		}
	};

	struct Point {
		double x;
		double y;
		double z;
		int32_t pointIndex;
		int32_t childIndex;
		double nx;
		double ny;
		double nz;
		float curvature;
	};

	// float CANCurvature(
	// 	const vector<Point>& points,
	// 	int64_t pointIndex,
	// 	double searchRadiusSq
	// ) {
	// 	auto point = points[pointIndex];
	// 	double px = point.x, py = point.y, pz = point.z,
	// 	npx = point.nx, npy = point.ny, npz = point.nz;
	// 	float curvature = 0.0f;
	// 	int neighbours = 0;
	//
	// 	for (int64_t j = 0; j < (int64_t)points.size(); j++) {
	// 		if (j == pointIndex) continue;
	//
	// 		double dx = points[j].x - px;
	// 		double dy = points[j].y - py;
	// 		double dz = points[j].z - pz;
	// 		double distSq = dx*dx + dy*dy + dz*dz;
	//
	// 		if (distSq > searchRadiusSq) continue;
	//
	// 		// project displacement onto tangent plane of p
	// 		double dot = dx * npx + dy * npy + dz * npz;
	// 		double tx = dx - dot * npx;
	// 		double ty = dy - dot * npy;
	// 		double tz = dz - dot * npz;
	//
	// 		double r2d = sqrt(tx*tx + ty*ty + tz*tz);
	// 		if (r2d < 1e-10) continue;
	//
	// 		float nqx = points[j].nx, nqy = points[j].ny, nqz = points[j].nz;
	//
	// 		float nxy = (float)((tx * nqx + ty * nqy + tz * nqz) / r2d);
	// 		float nz  = npx * nqx + npy * nqy + npz * nqz;
	//
	// 		float denom = sqrt(nxy * nxy + nz * nz) * (float)r2d;
	// 		if (denom < 1e-10f) continue;
	//
	// 		curvature += std::abs(-nxy / denom); //don't care about convex vs concave
	// 		neighbours++;
	// 	}
	//
	//    return neighbours > 0 ? curvature / neighbours : 0.0f;
	// }

	float gaussCurvature(
		const vector<Point>& points,
		int64_t pointIndex,
		double neighbourhoodRadiusSq
	) {
		auto point = points[pointIndex];
		double px = point.x, py = point.y, pz = point.z,
		npx = point.nx, npy = point.ny, npz = point.nz;

		vector<double> dist;
		dist.reserve(points.size());

		vector<size_t> neighbours;
		neighbours.reserve(maxNeighbours);


		HeapCmp cmp(dist);

		for (size_t i = 0; i < points.size(); ++i) {
			if (i == pointIndex) {
				dist.push_back(0);
				continue;
			}

			double dx = points[i].x - px;
			double dy = points[i].y - py;
			double dz = points[i].z - pz;
			double distSq = dx*dx + dy*dy + dz*dz;
			
			dist.push_back(distSq);
			
			if (distSq > neighbourhoodRadiusSq) 
				continue;
			
			if (neighbours.size() < maxNeighbours) {
				neighbours.push_back(i);
				push_heap(neighbours.begin(), neighbours.end(), cmp);
			} else {
				if (dist[neighbours.back()] > distSq) {
					neighbours.pop_back();
					pop_heap(neighbours.begin(), neighbours.end(), cmp);
					neighbours.push_back(i);
					push_heap(neighbours.begin(), neighbours.end(), cmp);
				}
			}
		}

		vector<Point> neighbourPoints;
		for (auto n : neighbours)
			neighbourPoints.push_back(points[n]);
		
		if (neighbourPoints.size() < 2) return 0.0f;

		thread_local std::mt19937 rng(std::random_device{}());
		std::uniform_int_distribution<size_t> distribution(0, neighbourPoints.size() - 1);

		double area = 0;
		double curvature = 0;

		for(size_t i = 0; i < samples; ++i) { // CNC-Uniform
			// pick 3 points at random from neighbourPoints
			size_t ai = distribution(rng);
			size_t bi = distribution(rng);
			size_t ci = distribution(rng);
			if (bi == ai || bi == ci || ai == ci) { --i; continue; }; // not a triangle, try again

			const Point& qa = neighbourPoints[ai];
			const Point& qb = neighbourPoints[bi];
			const Point& qc = neighbourPoints[ci];

			Vector3 xi = Vector3(qc.x, qc.y, qc.z);
			Vector3 ui = Vector3(qc.nx, qc.ny, qc.nz);

			Vector3 xj = Vector3(qa.x, qa.y, qa.z);
			Vector3 uj = Vector3(qa.nx, qa.ny, qa.nz);

			Vector3 xk = Vector3(qb.x, qb.y, qb.z);
			Vector3 uk = Vector3(qb.nx, qb.ny, qb.nz);

			auto uBar = (ui + uj + uk)/3;

			area += 0.5 * uBar.dot((xj - xi).cross(xk - xi));
			curvature += 0.5 * ui.dot(uj.cross(uk));
		}

		// if (std::abs(area) < 1e-10) return 0.0f;
		if (curvature < 0)
			curvature *= -1;

		return (float)(curvature / area);
	}

	// float meanCurvature(
	// 	const vector<Point>& points,
	// 	int64_t pointIndex,
	// 	double neighbourhoodRadiusSq
	// ) {
	// 	auto point = points[pointIndex];
	// 	double px = point.x, py = point.y, pz = point.z,
	// 	npx = point.nx, npy = point.ny, npz = point.nz;
	//
	// 	vector<double> dist;
	// 	dist.reserve(points.size());
	//
	// 	vector<size_t> neighbours;
	// 	neighbours.reserve(maxNeighbours);
	//
	// 	HeapCmp cmp(dist);
	//
	// 	for (size_t i = 0; i < points.size(); ++i) {
	// 		if (i == pointIndex) {
	// 			dist.push_back(0);
	// 			continue;
	// 		}
	//
	// 		double dx = points[i].x - px;
	// 		double dy = points[i].y - py;
	// 		double dz = points[i].z - pz;
	// 		double distSq = dx*dx + dy*dy + dz*dz;
	//
	// 		dist.push_back(distSq);
	//
	// 		if (distSq > neighbourhoodRadiusSq) 
	// 			continue;
	//
	// 		if (neighbours.size() < maxNeighbours) {
	// 			neighbours.push_back(i);
	// 			push_heap(neighbours.begin(), neighbours.end(), cmp);
	// 		} else {
	// 			if (dist[neighbours.back()] > distSq) {
	// 				neighbours.pop_back();
	// 				pop_heap(neighbours.begin(), neighbours.end(), cmp);
	// 				neighbours.push_back(i);
	// 				push_heap(neighbours.begin(), neighbours.end(), cmp);
	// 			}
	// 		}
	// 	}
	//
	// 	vector<Point> neighbourPoints;
	// 	for (auto n : neighbours)
	// 		neighbourPoints.push_back(points[n]);
	//
	// 	if (neighbourPoints.size() < 2) return 0.0f;
	//
	// 	thread_local std::mt19937 rng(std::random_device{}());
	// 	std::uniform_int_distribution<size_t> distribution(0, neighbourPoints.size() - 1);
	//
	// 	double area = 0;
	// 	double curvature = 0;
	//
	// 	for(size_t i = 0; i < samples; ++i) {
	// 		// pick 2 points at random from neighbourPoints
	// 		size_t ai = distribution(rng);
	// 		size_t bi = distribution(rng);
	// 		if (bi == ai) { --i; continue; };
	//
	// 		const Point& qj = neighbourPoints[ai];
	// 		const Point& qk = neighbourPoints[bi];
	//
	// 		Vector3 xi = Vector3(point.x, point.y, point.z);
	// 		Vector3 ui = Vector3(point.nx, point.ny, point.nz);
	//
	// 		Vector3 xj = Vector3(qj.x, qj.y, qj.z);
	// 		Vector3 uj = Vector3(qj.nx, qj.ny, qj.nz);
	//
	// 		Vector3 xk = Vector3(qk.x, qk.y, qk.z);
	// 		Vector3 uk = Vector3(qk.nx, qk.ny, qk.nz);
	//
	// 		auto uBar = (ui + uj + uk)/3;
	//
	// 		area += 0.5 * uBar.dot((xj - xi).cross(xk - xi));
	// 		// curvature += 0.5 * uBar.dot(
	// 		// 	(uk - uj).cross(xi)
	// 		// 	+ (ui - uk).cross(xj)
	// 		// 	+ (uj - ui).cross(xk)
	// 		// );
	// 		curvature += 0.5 * uBar.dot(
	// 			(uj - uk).cross(xj - xk) +
	// 			(uk - ui).cross(xk - xi) +
	// 			(ui - uj).cross(xi - xj)
	// 		);
	// 	}
	//
	// 	if (std::abs(area) < 1e-10) return 0.0f;
	// 	return (float)(curvature / area);
	//
	// 	return curvature / area;
	// }


	// subsample a local octree from bottom up
	void sample(Node* node, Attributes attributes, double baseSpacing, 
		function<void(Node*)> onNodeCompleted, 
		function<void(Node*)> onNodeDiscarded
	) {
		if (node->level() == 0) { // check for normals on root node
			bool hasNormals =
				attributes.get("NormalX") != nullptr &&
				attributes.get("NormalY") != nullptr &&
				attributes.get("NormalZ") != nullptr;

			if (!hasNormals) {
				cout << "ERROR: curvature sampler requires NormalX, NormalY and NormalZ.\n";
				cout << "Available attributes:\n";
				for (auto& a : attributes.list) {
					cout << "  " << a.name << "\n";
				}
				exit(1);
			}
		}
		int normalXOffset = attributes.getOffset("NormalX");
		int normalYOffset = attributes.getOffset("NormalY");
		int normalZOffset = attributes.getOffset("NormalZ");

		function<void(Node*, function<void(Node*)>)> traversePost = [&traversePost](Node* node, function<void(Node*)> callback) {
			for (auto child : node->children) {

				if (child != nullptr && !child->sampled) {
					traversePost(child.get(), callback);
				}
			}

			callback(node);
		};

		int64_t bytesPerPoint = attributes.bytes;
		Vector3 scale = attributes.posScale;
		Vector3 offset = attributes.posOffset;

		traversePost(node, [this, bytesPerPoint, baseSpacing, scale, offset, &onNodeCompleted, &onNodeDiscarded, attributes, normalXOffset, normalYOffset, normalZOffset](Node* node) {
			node->sampled = true;

			int64_t numPoints = node->numPoints;

			auto max = node->max;
			auto min = node->min;
			auto size = max - min;
			auto scale = attributes.posScale;
			auto offset = attributes.posOffset;

			bool isLeaf = node->isLeaf();

			if (isLeaf) {
				return false;
			}

			// =================================================================
			// SAMPLING
			// =================================================================
			//
			// first, check for each point whether it's accepted or rejected
			// save result in an array with one element for each point

			int64_t numPointsInChildren = 0;
			for (auto child : node->children) {
				if (child == nullptr) {
					continue;
				}

				numPointsInChildren += child->numPoints;
			}

			vector<Point> points;
			points.reserve(numPointsInChildren);

			vector<vector<int8_t>> acceptedChildPointFlags;
			vector<int64_t> numRejectedPerChild(8, 0);
			int64_t numAccepted = 0;

			for (int64_t childIndex = 0; childIndex < 8; childIndex++) {
				auto child = node->children[childIndex];

				if (child == nullptr) {
					acceptedChildPointFlags.push_back({});
					numRejectedPerChild.push_back({});

					continue;
				}

				vector<int8_t> acceptedFlags(child->numPoints, 0);
				acceptedChildPointFlags.push_back(acceptedFlags);

				for (int64_t i = 0; i < child->numPoints; i++) {
					int64_t pointOffset = i * attributes.bytes;
					int32_t* xyz = reinterpret_cast<int32_t*>(child->points->data_u8 + pointOffset);

					double x = (xyz[0] * scale.x) + offset.x;
					double y = (xyz[1] * scale.y) + offset.y;
					double z = (xyz[2] * scale.z) + offset.z;

					double nx = *reinterpret_cast<double*>(child->points->data_u8 + pointOffset + normalXOffset);
					double ny = *reinterpret_cast<double*>(child->points->data_u8 + pointOffset + normalYOffset);
					double nz = *reinterpret_cast<double*>(child->points->data_u8 + pointOffset + normalZOffset);

					Point point = { x, y, z, i, childIndex, nx, ny, nz, 0.0 };

					points.push_back(point);
				}

			}


			unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

			thread_local vector<Point> dbgAccepted(1'000'000);
			int64_t dbgNumAccepted = 0;
			double spacing = baseSpacing / pow(2.0, node->level());
			double squaredSpacing = spacing * spacing;
			double searchRadiusSq = squaredSpacing * 4.0;

			// calculate curvature
			for (int64_t i = 0; i < (int64_t)points.size(); i++) {
				auto& pi = points[i];
				pi.curvature = gaussCurvature(points, i, searchRadiusSq);
				std::cerr << 'C' << pi.curvature << '\n';
			}
			//normalize curvature
			float maxC = std::max_element(points.begin(), points.end(), 
								  [](const Point& a, const Point& b){ return a.curvature < b.curvature; })->curvature;
			if (maxC > 0) for (auto& p : points) p.curvature /= maxC;

			auto squaredDistance = [](Point& a, Point& b) {
				double dx = a.x - b.x;
				double dy = a.y - b.y;
				double dz = a.z - b.z;

				double dd = dx * dx + dy * dy + dz * dz;

				return dd;
			};

			auto center = (node->min + node->max) * 0.5;

			auto checkAccept = [/*&dbgChecks, &dbgSumChecks,*/ &dbgNumAccepted, spacing, squaredSpacing, &squaredDistance, center /*, &numDistanceChecks*/](Point candidate) {
				auto curvScale = candidate.curvature;
				auto curvSqSpacing = squaredSpacing * curvScale * curvScale;
				auto curvSpacing = spacing * curvScale;

				auto cx = candidate.x - center.x;
				auto cy = candidate.y - center.y;
				auto cz = candidate.z - center.z;
				auto cdd = cx * cx + cy * cy + cz * cz;
				auto cd = sqrt(cdd);
				auto limit = (cd - curvSpacing);
				auto limitSquared = limit * limit;

				int64_t j = 0;
				for (int64_t i = dbgNumAccepted - 1; i >= 0; i--) {

					auto& point = dbgAccepted[i];

					//dbgChecks++;
					//dbgSumChecks++;

					// check distance to center
					auto px = point.x - center.x;
					auto py = point.y - center.y;
					auto pz = point.z - center.z;
					auto pdd = px * px + py * py + pz * pz;
					//auto pd = sqrt(pdd);

					// stop when differences to center between candidate and accepted exceeds the spacing
					// any other previously accepted point will be even closer to the center.
					if (pdd < limitSquared) {
						return true;
					}

					double dd = squaredDistance(point, candidate);

					if (dd < curvSqSpacing) {
						return false;
					}

					j++;

					// also put a limit at x distance checks
					if (j > 10'000) {
						return true;
					}
				}

				return true;

			};

			auto parallel = std::execution::par_unseq;
			std::sort(parallel, points.begin(), points.end(), [center](Point a, Point b) -> bool {
				//sort by curvature first
				if (std::abs(a.curvature - b.curvature) > 1e-4f)
					return a.curvature > b.curvature; 

				auto ax = a.x - center.x;
				auto ay = a.y - center.y;
				auto az = a.z - center.z;
				auto add = ax * ax + ay * ay + az * az;

				auto bx = b.x - center.x;
				auto by = b.y - center.y;
				auto bz = b.z - center.z;
				auto bdd = bx * bx + by * by + bz * bz;

				// sort by distance to center
				return add < bdd;

				// sort by manhattan distance to center
				//return (ax + ay + az) < (bx + by + bz);

				// sort by z axis
				//return a.z < b.z;
			});

			for (Point point : points) {

				bool isAccepted = checkAccept(point);

				if (isAccepted) {
					dbgAccepted[dbgNumAccepted] = point;
					dbgNumAccepted++;
					numAccepted++;
				} else {
					numRejectedPerChild[point.childIndex]++;
				}

				acceptedChildPointFlags[point.childIndex][point.pointIndex] = isAccepted ? 1 : 0;


			}

			auto accepted = make_shared<Buffer>(numAccepted * attributes.bytes);
			for (int64_t childIndex = 0; childIndex < 8; childIndex++) {
				auto child = node->children[childIndex];

				if (child == nullptr) {
					continue;
				}

				auto numRejected = numRejectedPerChild[childIndex];
				auto& acceptedFlags = acceptedChildPointFlags[childIndex];
				auto rejected = make_shared<Buffer>(numRejected * attributes.bytes);

				for (int64_t i = 0; i < child->numPoints; i++) {
					auto isAccepted = acceptedFlags[i];
					int64_t pointOffset = i * attributes.bytes;

					if (isAccepted) {
						accepted->write(child->points->data_u8 + pointOffset, attributes.bytes);
						// rejected->write(child->points->data_u8 + pointOffset, attributes.bytes);
					} else {
						rejected->write(child->points->data_u8 + pointOffset, attributes.bytes);
					}
				}

				if (numRejected == 0 && child->isLeaf()) {
					onNodeDiscarded(child.get());

					node->children[childIndex] = nullptr;
				} if (numRejected > 0) {
					child->points = rejected;
					child->numPoints = numRejected;

					onNodeCompleted(child.get());
				} else if(numRejected == 0) {
					// the parent has taken all points from this child, 
					// so make this child an empty inner node.
					// Otherwise, the hierarchy file will claim that 
					// this node has points but because it doesn't have any,
					// decompressing the nonexistent point buffer fails
					// https://github.com/potree/potree/issues/1125
					child->points = nullptr;
					child->numPoints = 0;
					onNodeCompleted(child.get());
				}
			}

			node->points = accepted;
			node->numPoints = numAccepted;

			//{ // debug
			//	auto avgChecks = dbgSumChecks / points.size();
			//	string msg = "#checks: " + formatNumber(dbgSumChecks) + ", maxChecks: " + formatNumber(dbgMaxChecks) + ", avgChecks: " + formatNumber(avgChecks) + "\n";
			//	cout << msg;
			//}

			return true;
		});
	}

};

