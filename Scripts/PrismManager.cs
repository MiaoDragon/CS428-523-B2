using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;




public class PrismManager : MonoBehaviour
{
    public int prismCount = 10;
    public float prismRegionRadiusXZ = 5;
    public float prismRegionRadiusY = 5;
    public float maxPrismScaleXZ = 5;
    public float maxPrismScaleY = 5;
    public GameObject regularPrismPrefab;
    public GameObject irregularPrismPrefab;

    private List<Prism> prisms = new List<Prism>();
    private List<GameObject> prismObjects = new List<GameObject>();
    private GameObject prismParent;
    private Dictionary<Prism,bool> prismColliding = new Dictionary<Prism, bool>();

    private const float UPDATE_RATE = 0.5f;

    #region Unity Functions

    void Start()
    {
        Random.InitState(0);    //10 for no collision

        prismParent = GameObject.Find("Prisms");
        for (int i = 0; i < prismCount; i++)
        {
            var randPointCount = Mathf.RoundToInt(3 + Random.value * 7);
            var randYRot = Random.value * 360;
            var randScale = new Vector3((Random.value - 0.5f) * 2 * maxPrismScaleXZ, (Random.value - 0.5f) * 2 * maxPrismScaleY, (Random.value - 0.5f) * 2 * maxPrismScaleXZ);
            var randPos = new Vector3((Random.value - 0.5f) * 2 * prismRegionRadiusXZ, (Random.value - 0.5f) * 2 * prismRegionRadiusY, (Random.value - 0.5f) * 2 * prismRegionRadiusXZ);

            GameObject prism = null;
            Prism prismScript = null;
            if (Random.value < 0.5f)
            {
                prism = Instantiate(regularPrismPrefab, randPos, Quaternion.Euler(0, randYRot, 0));
                prismScript = prism.GetComponent<RegularPrism>();
            }
            else
            {
                prism = Instantiate(irregularPrismPrefab, randPos, Quaternion.Euler(0, randYRot, 0));
                prismScript = prism.GetComponent<IrregularPrism>();
            }
            prism.name = "Prism " + i;
            prism.transform.localScale = randScale;
            prism.transform.parent = prismParent.transform;
            prismScript.pointCount = randPointCount;
            prismScript.prismObject = prism;

            prisms.Add(prismScript);
            prismObjects.Add(prism);
            prismColliding.Add(prismScript, false);
        }

        StartCoroutine(Run());
    }
    
    void Update()
    {
        #region Visualization

        DrawPrismRegion();
        DrawPrismWireFrames();

#if UNITY_EDITOR
        if (Application.isFocused)
        {
            UnityEditor.SceneView.FocusWindowIfItsOpen(typeof(UnityEditor.SceneView));
        }
#endif

        #endregion
    }

    IEnumerator Run()
    {
        yield return null;

        while (true)
        {
            foreach (var prism in prisms)
            {
                prismColliding[prism] = false;
            }

            foreach (var collision in PotentialCollisions())
            {
                if (CheckCollision(collision))
                {
                    prismColliding[collision.a] = true;
                    prismColliding[collision.b] = true;

                    ResolveCollision(collision);
                }
            }

            yield return new WaitForSeconds(UPDATE_RATE);
        }
    }

    #endregion

    #region Incomplete Functions


    private class BoundingBoxAxis
    {
        public float value;
        public int object_id;
        public int type;
        public BoundingBoxAxis(float val, int oid, int tp)
        {
            value = val;
            object_id = oid;
            type = tp;
        }        
    }

    private IEnumerable<PrismCollision> PotentialCollisions()
    {
        /*
         * default implementation
        for (int i = 0; i < prisms.Count; i++) {
            for (int j = i + 1; j < prisms.Count; j++) {
                var checkPrisms = new PrismCollision();
                checkPrisms.a = prisms[i];
                checkPrisms.b = prisms[j];

                yield return checkPrisms;
            }
        }

        yield break;
        */

        // * implementation using sort & sweep
        // for each prism shape, obtain the bounding box
        // by checking each point, and obtain the min and max
        // in each dimension


        // *TODO: the 3D points are not the points on the polygons, but something different
        // need to change the algorithm to take care of that

        bool[] potentialObjects = new bool[prisms.Count];
        // counting how many times the pair appeared
        int[] potentialPairHashes = new int[prisms.Count * prisms.Count];
        // hashing function: big_idx*N + small_idx

        List<Vector3[]> prismPoints = new List<Vector3[]>();
        for (int i=0; i< prisms.Count(); i++)
        {
            if (maxPrismScaleY != 0)
            {
                Vector3[] points = Obtain3DPoints(prisms[i]);
                prismPoints.Add(points);
            }
            else
            {
                prismPoints.Add(prisms[i].points);
            }
        }
        // initialize the potential object list with all the objects in the scene
        for (int i = 0; i < prisms.Count; i++) potentialObjects[i] = true;
        for (int i = 0; i < prisms.Count * prisms.Count; i++) potentialPairHashes[i] = 0;
        for (int dim=0; dim<3; dim++)
        {
            //Debug.Log("dim: " + dim.ToString());
            // skip the second dimension if the object is 2D
            if (maxPrismScaleY==0 && dim==1)
            {
                //Debug.Log("skipping");

                // one of the dimension is 0. Add count 1 to each of the potential pair
                for (int i = 0; i < prisms.Count * prisms.Count; i++) potentialPairHashes[i] += 1;
                continue;
            }
            List<BoundingBoxAxis> v_list = new List<BoundingBoxAxis>(prisms.Count);
            // value, obj_idx, {0,1} indicating start/end
            for (int i = 0; i < prisms.Count; i++)
            {
                if (!potentialObjects[i]) continue;
                float v_min = Mathf.Max(prismRegionRadiusXZ+maxPrismScaleXZ*2, prismRegionRadiusY+maxPrismScaleY*2);
                float v_max = -v_min;
                for (int j = 0; j < prismPoints[i].Length; j++)
                {
                    if (prismPoints[i][j][dim] < v_min)
                    {
                        v_min = prismPoints[i][j][dim];
                    }
                    if (prismPoints[i][j][dim] > v_max)
                    {
                        v_max = prismPoints[i][j][dim];
                    }
                }
                //Debug.Log("object " + i.ToString() + " vmin: " + v_min.ToString());
                //Debug.Log("object " + i.ToString() + " vmax: " + v_max.ToString());

                // add corresponding axis to the array
                v_list.Add(new BoundingBoxAxis(v_min, i, 0));
                v_list.Add(new BoundingBoxAxis(v_max, i, 1));
            }
            // sort the list by the x/y/z values

            v_list.Sort((p1, p2) => p1.value.CompareTo(p2.value));

            //Debug.Log("after sorting...");

            List<int> activeObjects = new List<int>(v_list.Count());

            // sweep to obtain the pair
            // reset potential objects
            for (int i=0; i < prisms.Count; i++) potentialObjects[i] = false;
            for (int i=0; i<v_list.Count(); i++)
            {
                // if i-th element is the start, then check it with existing active objects
                //      then add it to the active set
                int oid = v_list[i].object_id;
                int type = v_list[i].type;
                if (type == 0) 
                {
                    for (int j = 0; j < activeObjects.Count(); j++)
                    {
                        // the two objects are potential pairs, check if we encountered them
                        // before for "dim" times. If not, then they're not potential pairs
                        // use hashing to check if the pair appeared before
                        int smaller = Mathf.Min(activeObjects[j], oid);
                        int bigger = Mathf.Max(activeObjects[j], oid);

                        // if the pair didn't appear enough times before, they'll be ruled out
                        if (!(potentialPairHashes[bigger * prisms.Count + smaller] == dim)) continue;
                        // potential pair will add one to count
                        potentialPairHashes[bigger * prisms.Count + smaller] += 1;
                        
                        potentialObjects[activeObjects[j]] = true;
                        potentialObjects[oid] = true;

                        // yield the current pair if the pair has appeared for enough times (3 times for dim 3)
                        if (potentialPairHashes[bigger * prisms.Count + smaller] == 3)
                        {
                            // yield the current pair
                            var checkPrisms = new PrismCollision();
                            checkPrisms.a = prisms[oid];
                            checkPrisms.b = prisms[activeObjects[j]];
                            //Debug.Log("collision pairs: object " + oid.ToString() + ", object " + activeObjects[j].ToString());
                            yield return checkPrisms;
                        }
                    }
                    activeObjects.Add(oid);
                }
                // if i-th element is the end, then remove it from the active object set
                else
                {
                    activeObjects.Remove(oid);
                }
            }
        }
        yield break;
    }


    private Vector3[] MinkowskiDiff(Vector3[] pointsA, Vector3[] pointsB)
    {
        // obtain the MinkowskiDiff of two convex shapes
        // since SupportPoint will always return the point no the edge of the convex shape, we can 
        // just return the point set after MinkowskiDiff
        Vector3[] points = new Vector3[pointsA.Length*pointsB.Length];
        for (int i = 0; i < pointsA.Length; i++)
        {
            for (int j = 0; j < pointsB.Length; j++)
            {
                points[i * pointsB.Length + j] = pointsA[i] - pointsB[j];
            }
        }

        return points;
    }

    private int GetSupportPointIdx(Vector3 v, Vector3[] points, Vector3 origin)
    {
        // get the support point along vector (origin-v) relative to origin from points
        // output the index of the point
        float minValue = (maxPrismScaleXZ + prismRegionRadiusXZ) * 5;
        int minIdx = -1;
        Vector3 direction = (v - origin) / (v - origin).magnitude;

        for (int i=0; i< points.Length; i++)
        {
            float distance = Vector3.Dot(direction, points[i] - origin);
            if (distance < minValue)
            {
                minValue = distance;
                minIdx = i;
            }
        }
        return minIdx;
    }

    private float LinePointDistance(Vector3 a, Vector3 b, Vector3 p)
    {
        // compute the distance from p to line a-b
        Vector3 ab = b - a;
        // assume we are in the same Y plane
        Vector3 normal = new Vector3(-ab[2], 0, ab[0]);
        normal = normal / normal.magnitude;
        float distance = Vector3.Dot(p - a, normal);
        return Mathf.Abs(distance);
    }

    private Vector3 LineNormal(Vector3 a, Vector3 b, Vector3 c)
    {
        // obtain the line normal of ab, with c being the third point
        Vector3 ab = b - a;
        Vector3 ac = c - a;
        Vector3 normal_ab = new Vector3(-ab[2], 0, ab[0]);
        normal_ab = normal_ab / normal_ab.magnitude;
        if (Vector3.Dot(ac, normal_ab) > 0)
        {
            normal_ab = -normal_ab;
        }
        return normal_ab;
    }

    private bool TriangleInside(Vector3 a, Vector3 b, Vector3 c, Vector3 p)
    {
        // obtain the normal for each line
        Vector3 normal_ab = LineNormal(a, b, c);
        Vector3 normal_ac = LineNormal(a, c, b);
        Vector3 normal_bc = LineNormal(b, c, a);
        // check if p lies in the triangle
        Vector3 ap = p - a;
        Vector3 bp = p - b;



        if (Vector3.Dot(normal_ab, ap) <= 1e-8f && Vector3.Dot(normal_ac, ap) <= 1e-8f && Vector3.Dot(normal_bc, bp) <= 1e-8f)
        {
            return true;
        }
        return false;
    }

    private Vector3 ClosestPointToSimplex(List<Vector3> simplex, Vector3 p)
    {
        // dimension only up to 3 (triangle)
        // for each dimension, see if it is the min distance
        Vector3 minPoint = simplex[0];
        float minValue = (minPoint - p).magnitude;
        // dimension 1: point
        for (int i=0; i<simplex.Count; i++)
        {
            if ((simplex[i] - p).magnitude < minValue)
            {
                minPoint = simplex[i];
                minValue = (simplex[i] - p).magnitude;
            }
        }
        for (int i=0; i<simplex.Count; i++)
        {
            for (int j=i+1; j<simplex.Count; j++)
            {
                // if the projection is inside the line
                Vector3 a = simplex[i];
                Vector3 b = simplex[j];
                if (Vector3.Dot(b-a, p-a)>0 && Vector3.Dot(a - b, p - b) > 0)
                {
                    Vector3 direction = (b - a) / (b - a).magnitude;
                    Vector3 project = Vector3.Dot(direction, p - a)*direction + a;
                    if ((project - p).magnitude < minValue)
                    {
                        minPoint = project;
                        minValue = (project - p).magnitude;
                    }
                }
            }
        }
        if (simplex.Count == 3)
        {
            // final case: check with triangle
            // for triangle we can form an equation and solve it.
            // assume the perpendicular point is at: alpha * (c-a) + beta * (b-a)
            // then we have:
            // (alpha * (c-a) + beta * (b-a) - p)T(c-a) = 0
            // (alpha * (c-a) + beta * (b-a) - p)T(b-a) = 0
            Vector3 a = simplex[0];
            Vector3 b = simplex[1];
            Vector3 c = simplex[2];
            Vector3 ab = b - a;
            Vector3 ac = c - a;
            float alpha = 0.0f, beta = 0.0f;
            if (Vector3.Dot(ab, ac) == 0)
            {
                // degenerate case
                alpha = Vector3.Dot(p, ac) / Vector3.Dot(ac, ac);
                beta = Vector3.Dot(p, ab) / Vector3.Dot(ab, ab);                
            }
            else
            {
                beta = (Vector3.Dot(p, ab) * Vector3.Dot(ac, ac) - Vector3.Dot(p, ac) * Vector3.Dot(ac, ab)) /
                        (Vector3.Dot(ab, ab) * Vector3.Dot(ac, ac) - Vector3.Dot(ab, ac) * Vector3.Dot(ab, ac));
                alpha = (Vector3.Dot(p, ac) * Vector3.Dot(ab, ab) - Vector3.Dot(p, ab) * Vector3.Dot(ac, ab)) /
                        (Vector3.Dot(ac, ac) * Vector3.Dot(ab, ab) - Vector3.Dot(ab, ac) * Vector3.Dot(ab, ac));
            }
            if (beta >= 0.0f && beta <= 1.0f && alpha >= 0.0f && alpha <= 1.0f)
            {
                Vector3 project = alpha * ac + beta * ab;
                if ((project - p).magnitude < minValue)
                {
                    minPoint = project;
                    minValue = (project - p).magnitude;
                }
            }

        }
        return minPoint;

    }


    private void ClosestPointToEdges(List<Vector3> edgeA, List<Vector3> edgeB, Vector3 p, 
                                                out Vector3 v, out int edgeIdx)
    {
        /* return the edgeA, minPoint, edgeB */
        // dimension only up to 3 (triangle)
        // for each dimension, see if it is the min distance
        Vector3 minPoint = edgeA[0];
        int minIdx = 0;
        float minValue = (minPoint - p).magnitude;
        // dimension 1: point
        for (int i = 0; i < edgeA.Count; i++)
        {
            if ((edgeA[i] - p).magnitude < minValue)
            {
                minPoint = edgeA[i];
                minValue = (edgeA[i] - p).magnitude;
                minIdx = i;
            }
        }

        for (int i = 0; i < edgeB.Count; i++)
        {
            if ((edgeB[i] - p).magnitude < minValue)
            {
                minPoint = edgeB[i];
                minValue = (edgeB[i] - p).magnitude;
                minIdx = i;
            }
        }

        for (int edge_i = 0; edge_i < edgeA.Count; edge_i++)
        {
            // if the projection is inside the line
            Vector3 a = edgeA[edge_i];
            Vector3 b = edgeB[edge_i];
            if (Vector3.Dot(b - a, p - a) > 0 && Vector3.Dot(a - b, p - b) > 0)
            {
                Vector3 direction = (b - a) / (b - a).magnitude;
                Vector3 project = Vector3.Dot(direction, p - a) * direction + a;
                if ((project - p).magnitude < minValue)
                {
                    minPoint = project;
                    minValue = (project - p).magnitude;
                    minIdx = edge_i;
                }
            }
        }
        v = minPoint;
        edgeIdx = minIdx;
    }


    private Vector3 EPA2D(List<Vector3> W, Vector3[] points)
    {
        Vector3 origin = new Vector3(0.0f, 0.0f, 0.0f);
        List<Vector3> edgeA = new List<Vector3>();
        List<Vector3> edgeB = new List<Vector3>();
        for (int i=0; i<W.Count(); i++)
        {
            for (int j=i+1; j<W.Count(); j++)
            {
                edgeA.Add(W[i]);
                edgeB.Add(W[j]);
            }
        }
        Vector3 supportPoint = points[0];
        while (true)
        {
            // find closest point of simplex to origin
            Vector3 v = new Vector3();
            int edgeIdx = 0;
            ClosestPointToEdges(edgeA, edgeB, origin, out v, out edgeIdx);
            // find the supporting point in the Minkowski point set for v

            int supportIdx = GetSupportPointIdx(2 * origin - v, points, origin);
            Vector3 newSupportPoint = points[supportIdx];

            if (Vector3.Distance(supportPoint, newSupportPoint) <= 1e-4f)
            {
                return supportPoint;
            }
            supportPoint = newSupportPoint;
            Vector3 a = edgeA[edgeIdx];
            Vector3 b = edgeB[edgeIdx];
            edgeA.RemoveAt(edgeIdx);
            edgeB.RemoveAt(edgeIdx);
            edgeA.Add(a);
            edgeB.Add(supportPoint);
            edgeA.Add(supportPoint);
            edgeB.Add(b);
        }

    }

    private void GJKEPA2D(Vector3[] points, out bool inCollision, out Vector3 penetrationDepth)
    {
        /* return -1 if not incollision. otherwise return the positive penetration depth*/
   
        // * GJK algorithm for collision checking
        Vector3 v = points[0];
        Vector3 origin = new Vector3(0.0f, 0.0f, 0.0f);
        List<Vector3> W = new List<Vector3>();
        W.Add(v);
        int secondIdx = GetSupportPointIdx(v, points, origin);
        W.Add(points[secondIdx]);
        inCollision = false;
        float threshold = 1e-4f;
        while (true)
        {
            // find closest point from the simplex to origin
            Vector3 v_new = ClosestPointToSimplex(W, origin);

            // check convergence
            if (Mathf.Abs(v_new.magnitude - v.magnitude) <= threshold)
            {
                inCollision = false;
                break;
            }
            // find support point based on the new vector
            int supportIdx = GetSupportPointIdx(v_new, points, origin);

            // put the new point to W
            //Debug.Log("v_new: " + v_new.ToString());
            W.Add(points[supportIdx]);

            // check if the origin is in the convex hull of W
            // check by seeing if the point lies inside the line for each of the line
            if (W.Count() < 3 && Vector3.Distance(ClosestPointToSimplex(W, origin), origin) <= 1e-8f)
            {
                // edge case: when the origin is on the edge
                inCollision = true;
                break;
            }
            if (W.Count() == 3 && TriangleInside(W[0], W[1], W[2], origin))
            {
                // converged, and inside
                inCollision = true;
                break;
            }

            // remove points from W that are not the closest
            // for 2D check each line, and find the two points that correspond to the min-distance
            // remove the other point

            if (W.Count() == 3)
            {
                float distance1 = LinePointDistance(W[0], W[1], origin);
                float distance2 = LinePointDistance(W[0], W[2], origin);
                float distance3 = LinePointDistance(W[1], W[2], origin);
                float minDistance = Mathf.Min(distance1, distance2, distance3);
                if (minDistance == distance1) W.RemoveAt(2);
                else if (minDistance == distance2) W.RemoveAt(1);
                else W.RemoveAt(0);
            }
            v = v_new;
        }

        // * after GJK algorithm, obtain penetration depth from EPA
        if (!inCollision) penetrationDepth = Vector3.zero;
        else
        {
            // EPA algorithm for finding the penetration depth
            penetrationDepth = EPA2D(W, points);
        }
    }

    private bool boundaryCheck2D(Vector3[] points, Vector3 origin)
    {
        // check whether the origin lies on the boundary of ConvexHull(points)
        // just need to check if the edge connecting origin and one of the points
        // can support all the points on one side
        for (int i=0; i<points.Length; i++)
        {
            if (Vector3.Distance(points[i], origin) <= 1e-7f)
            {
                continue;
                // close to zero
            }
            Vector3 edge = points[i] - origin;
            edge = edge / edge.magnitude;
            Vector3 normal = new Vector3(-edge[2], edge[1], edge[0]);
            int numNegative = 0;
            int numPositive = 0;
            int numEqual = 0;
            bool supported = true;
            for (int j=0; j<points.Length; j++)
            {
                if (Vector3.Dot(points[j] - origin, normal) > 1e-7f)
                {
                    numPositive += 1;
                }
                else if (Vector3.Dot(points[j] - origin, normal) > -1e-7f)
                {
                    numEqual += 1;
                }
                else
                {
                    numNegative += 1;
                }
                if (numPositive > 0 && numNegative > 0)
                {
                    supported = false;
                    break;
                }
            }
            if (supported)
            {
                return true;
            }
        }
        return false;
    }

    private bool CheckCollision2D(PrismCollision collision)
    {
        var prismA = collision.a;
        var prismB = collision.b;

        collision.penetrationDepthVectorAB = Vector3.zero;

        Vector3[] points = MinkowskiDiff(prismA.points, prismB.points);

        // boundary check if the origin lies on the boundary of the minkowskiDiff
        if (boundaryCheck2D(points, Vector3.zero))
        {
            return false;  // in contact with each other, but with zero penetration vector
        }

        // otherwise the origin lies in the interior or exterior
        Vector3 penetrationDepth = new Vector3();
        bool inCollision = false;
        GJKEPA2D(points, out inCollision, out penetrationDepth);
        if (!inCollision)
        {
            return false;
        }
        else
        {
            collision.penetrationDepthVectorAB = penetrationDepth;
        }

        return inCollision;

    }

    private Vector3[] Obtain3DPoints(Prism prism)
    {
        // obtain the 3D points (top and bottom) from prism
        var yMin = prism.midY - prism.height / 2 * prism.prismObject.transform.localScale.y;
        var yMax = prism.midY + prism.height / 2 * prism.prismObject.transform.localScale.y;
        Vector3[] points = new Vector3[prism.points.Length * 2];
        for (int i=0; i < prism.points.Length; i++)
        {
            points[2*i] = prism.points[i] + Vector3.up * yMin;
            points[2 * i + 1] = prism.points[i] + Vector3.up * yMax;
        }
        return points;
    }
    private bool CheckCollision3D(PrismCollision collision)
    {
        // For 3D, we firstly check if the minkowski points in the z axis contain origin
        // if so, we just need to do minkowskidiff for one plane, and then check
        // the Y-plane where Y=0
        // The penetration vector will also be horizontal since this is the shortest vector
        var prismA = collision.a;
        var prismB = collision.b;

        collision.penetrationDepthVectorAB = Vector3.zero;

        // obtain pointsA and pointsB from the prism
        Vector3[] pointsA = Obtain3DPoints(prismA);
        Vector3[] pointsB = Obtain3DPoints(prismB);

        Vector3[] points = MinkowskiDiff(pointsA, pointsB);

        // pre-check: if the max_y >0 and min_y < 0
        float maxY = -1.0f;
        int maxYCount = 0;
        float minY = 1.0f;
        for (int i=0; i<points.Length; i++)
        {
            if (points[i].y > maxY)
            {
                maxY = points[i].y;
                maxYCount = 0;
            }
            if (points[i].y < minY)
            {
                minY = points[i].y;
            }
            if (points[i].y == maxY)
            {
                maxYCount += 1;
            }
        }
        //Debug.Log("maxY: " + maxY.ToString());
        //Debug.Log("minY: " + minY.ToString());

        if (maxY < 0 || minY > 0)
        {
            // not in collision
            return false;
        }

        // otherwise, obtain the projection on Y-plane: Y=0
        Vector3[] newPoints = new Vector3[maxYCount];
        int count = 0;
        for (int i=0; i<points.Length; i++)
        {
            if (points[i].y == maxY)
            {
                newPoints[count].x = points[i].x;
                newPoints[count].y = 0;
                newPoints[count].z = points[i].z;
                count += 1;
            }
        }

        // boundary check if the origin lies on the boundary of the minkowskiDiff
        if (boundaryCheck2D(newPoints, Vector3.zero))
        {
            return false;  // in contact with each other, but with zero penetration vector
        }

        // otherwise the origin lies in the interior or exterior
        Vector3 penetrationDepth = new Vector3();
        bool inCollision = false;
        GJKEPA2D(newPoints, out inCollision, out penetrationDepth);
        if (!inCollision)
        {
            return false;
        }
        else
        {
            collision.penetrationDepthVectorAB = penetrationDepth;
        }

        return inCollision;

    }

    private bool CheckCollision(PrismCollision collision)
    {
        var prismA = collision.a;
        var prismB = collision.b;
        if (maxPrismScaleY > 0)
            return CheckCollision3D(collision);
        else
            return CheckCollision2D(collision);

        collision.penetrationDepthVectorAB = Vector3.zero;

        return true;

    }

    #endregion

    #region Private Functions

    private void ResolveCollision(PrismCollision collision)
    {
        var prismObjA = collision.a.prismObject;
        var prismObjB = collision.b.prismObject;

        var pushA = -collision.penetrationDepthVectorAB / 2;
        var pushB = collision.penetrationDepthVectorAB / 2;

        for (int i = 0; i < collision.a.pointCount; i++)
        {
            collision.a.points[i] += pushA;
        }
        for (int i = 0; i < collision.b.pointCount; i++)
        {
            collision.b.points[i] += pushB;
        }
        //prismObjA.transform.position += pushA;
        //prismObjB.transform.position += pushB;

        Debug.DrawLine(prismObjA.transform.position, prismObjA.transform.position + collision.penetrationDepthVectorAB, Color.cyan, UPDATE_RATE);
    }
    
    #endregion

    #region Visualization Functions

    private void DrawPrismRegion()
    {
        var points = new Vector3[] { new Vector3(1, 0, 1), new Vector3(1, 0, -1), new Vector3(-1, 0, -1), new Vector3(-1, 0, 1) }.Select(p => p * prismRegionRadiusXZ).ToArray();
        
        var yMin = -prismRegionRadiusY;
        var yMax = prismRegionRadiusY;

        var wireFrameColor = Color.yellow;

        foreach (var point in points)
        {
            Debug.DrawLine(point + Vector3.up * yMin, point + Vector3.up * yMax, wireFrameColor);
        }

        for (int i = 0; i < points.Length; i++)
        {
            Debug.DrawLine(points[i] + Vector3.up * yMin, points[(i + 1) % points.Length] + Vector3.up * yMin, wireFrameColor);
            Debug.DrawLine(points[i] + Vector3.up * yMax, points[(i + 1) % points.Length] + Vector3.up * yMax, wireFrameColor);
        }
    }

    private void DrawPrismWireFrames()
    {
        for (int prismIndex = 0; prismIndex < prisms.Count; prismIndex++)
        {
            var prism = prisms[prismIndex];
            var prismTransform = prismObjects[prismIndex].transform;

            var yMin = prism.midY - prism.height / 2 * prismTransform.localScale.y;
            var yMax = prism.midY + prism.height / 2 * prismTransform.localScale.y;

            var wireFrameColor = prismColliding[prisms[prismIndex]] ? Color.red : Color.green;

            foreach (var point in prism.points)
            {
                Debug.DrawLine(point + Vector3.up * yMin, point + Vector3.up * yMax, wireFrameColor);
            }

            for (int i = 0; i < prism.pointCount; i++)
            {
                Debug.DrawLine(prism.points[i] + Vector3.up * yMin, prism.points[(i + 1) % prism.pointCount] + Vector3.up * yMin, wireFrameColor);
                Debug.DrawLine(prism.points[i] + Vector3.up * yMax, prism.points[(i + 1) % prism.pointCount] + Vector3.up * yMax, wireFrameColor);
            }
        }
    }

    #endregion

    #region Utility Classes

    private class PrismCollision
    {
        public Prism a;
        public Prism b;
        public Vector3 penetrationDepthVectorAB;
    }

    private class Tuple<K,V>
    {
        public K Item1;
        public V Item2;

        public Tuple(K k, V v) {
            Item1 = k;
            Item2 = v;
        }
    }

    #endregion
}
