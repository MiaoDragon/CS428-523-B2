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
        bool[] potentialObjects = new bool[prisms.Count];
        // counting how many times the pair appeared
        int[] potentialPairHashes = new int[prisms.Count * prisms.Count];
        // hashing function: big_idx*N + small_idx

        // initialize the potential object list with all the objects in the scene
        for (int i = 0; i < prisms.Count; i++) potentialObjects[i] = true;
        for (int i = 0; i < prisms.Count * prisms.Count; i++) potentialPairHashes[i] = 0;
        for (int dim=0; dim<3; dim++)
        {
            // skip the second dimension if the object is 2D
            if (maxPrismScaleY==0 && dim==1)
            {
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
                for (int j = 0; j < prisms[i].pointCount; j++)
                {
                    if (prisms[i].points[j][dim] < v_min)
                    {
                        v_min = prisms[i].points[j][dim];
                    }
                    if (prisms[i].points[j][dim] > v_max)
                    {
                        v_max = prisms[i].points[j][dim];
                    }
                }
                // add corresponding axis to the array
                v_list.Add(new BoundingBoxAxis(v_min, i, 0));
                v_list.Add(new BoundingBoxAxis(v_max, i, 1));
            }
            // sort the list by the x/y/z values

            v_list.Sort((p1, p2) => p1.value.CompareTo(p2.value));

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

    private bool CheckCollision(PrismCollision collision)
    {
        var prismA = collision.a;
        var prismB = collision.b;

        
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
