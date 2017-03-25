// Define all the PhotoMapCollection routines.

#include "PhotoMapCollection.h"
#include "StringStuff.h"

using namespace photometry;

// The word that we expect all serialized collections to have on their first line
const string
PhotoMapCollection::magicKey = "PhotoMapCollection";

//////////////////////////////////////////////////////////////
// Static data
//////////////////////////////////////////////////////////////
map<string,PhotoMapCollection::Creator>
PhotoMapCollection::mapCreators;

bool
PhotoMapCollection::creatorsInitialized=false;

void
PhotoMapCollection::PhotoMapTypeInitialize() {
  if (creatorsInitialized) return;
  creatorsInitialized = true;
  registerMapType<IdentityMap>();
  registerMapType<ConstantMap>();
  registerMapType<PolyMap>();
}

//////////////////////////////////////////////////////////////
// Construction / destruction
//////////////////////////////////////////////////////////////

PhotoMapCollection::PhotoMapCollection() {
  PhotoMapTypeInitialize();
  rebuildParameterVector();
}

PhotoMapCollection::~PhotoMapCollection() {
  // Destroy all of the compounded maps and submaps we've made:
  for (MapIter i = mapElements.begin(); i!=mapElements.end(); ++i) {
    if (i->second.atom) {
      delete i->second.atom;
      i->second.atom = 0;
    }
    if (i->second.realization) {
      delete i->second.realization;
      i->second.realization = 0;
    }
  }
}

//////////////////////////////////////////////////////////////
// Parameter manipluation
//////////////////////////////////////////////////////////////

void
PhotoMapCollection::setParams(const DVector& p) {
  Assert(p.size()==parameterCount);
  for (MapIter i = mapElements.begin(); i!=mapElements.end(); ++i) {
    MapElement& map = i->second;
    if (map.isFixed) continue;
    int nSub = map.nParams;
    if (nSub<=0) continue;
    DVector subp(nSub, 0.);
    subp = p.subVector(map.startIndex, map.startIndex+nSub);
    Assert(map.atom)
    map.atom->setParams(subp);
  }
}

DVector
PhotoMapCollection::getParams() const {
  DVector p(parameterCount, 888.);
  for (auto& melpair : mapElements) {
    const MapElement& map = melpair.second;
    if (map.isFixed) continue;
    int nSub = map.nParams;
    if (nSub<=0) continue;
    if (!map.atom) 
      cerr << "mapElement is not atomic: " << melpair.first << " params: " << nSub << endl;
    Assert(map.atom);
    p.subVector(map.startIndex, map.startIndex+nSub) = 
      map.atom->getParams().subVector(0,nSub);
  }
  return p;
}

void
PhotoMapCollection::copyParamsFrom(const PhotoMap& pm) {
  auto melpair = mapElements.find(pm.getName());
  if (melpair == mapElements.end())
    throw PhotometryError("Attempt to copyParamsFrom non-existent map element " + pm.getName());
  if (!melpair->second.atom)
    throw PhotometryError("Attempt to copyParamsFrom non-atomic map element " + pm.getName());
  melpair->second.atom->setParams(pm.getParams());
}

void 
PhotoMapCollection::setFixed(string name, bool isFixed) {
  set<string> s;
  s.insert(name);
  setFixed(s, isFixed);
}

void 
PhotoMapCollection::setFixed(set<string> nameList, bool isFixed) {
  set<string> fixThese;
  for (auto iName : nameList) {
    set<string> addThese = dependencies(iName);
    fixThese.insert(addThese.begin(), addThese.end());
  }
  for (auto iName : fixThese) {
    auto& mel = mapElements[iName];
    // Update isFixed if this atom has free parameters
    if (mel.atom && mel.atom->nParams()>0) {
      mel.isFixed = isFixed;
    }
  }
  // Now calculate new parameter indices and propagate new configuration to all SubMaps
  rebuildParameterVector();
}

bool
PhotoMapCollection::getFixed(string name) const {
  ConstMapIter i = mapElements.find(name);
  if (i==mapElements.end())
    throw PhotometryError("No member of PhotoMapCollection named " + name);
  auto& mel = i->second;
  if (mel.atom) {
    return mel.isFixed;
  } else {
    for (auto i : mel.subordinateMaps) {
      if (!getFixed(i)) {
	// Any free member means compound is free
	return false;
      }
    }
    return true;  // Fixed if all subordinates are
  }
}

bool
PhotoMapCollection::isAtomic(string name) const {
  auto i = mapElements.find(name);
  if (i==mapElements.end())
    throw PhotometryError("No member of PhotoMapCollection named " + name);
  if (i->second.atom) {
    return true;
  } else {
    return false;
  }
}

void
PhotoMapCollection::rebuildParameterVector() {
  // First assign new starting indices to every atomic map element

  // Restart the parameter counting:
  parameterCount = 0;
  // And map counting
  atomCount = 0;
  freeCount = 0;
  for (MapIter i = mapElements.begin(); i!=mapElements.end(); ++i) {
    MapElement& map = i->second;
    // Only atomic map components go into the big parameter vector
    map.nParams = 0;
    map.number = -1;
    if (!(map.atom)) continue;
    atomCount++;
    int nMap = map.atom->nParams();
    if (nMap>0 && !map.isFixed) {
      // Have some parameters; append them to master list.
      map.nParams = nMap;
      map.startIndex = parameterCount;
      // Each atom with free parameters gets a serial number.
      map.number = freeCount;
      parameterCount += nMap;
      ++freeCount;
    }
  }

  // Then go through all extant SubMaps and update their versions of vectors & parameter counts.
  for (MapIter i = mapElements.begin(); i!=mapElements.end(); ++i) {
    MapElement& map = i->second;
    // Is there a SubMap here that we need to update?
    SubMap* sub = map.realization;
    if (!sub) continue;

    // Loop through all the map elements that this SubMap uses
    for (int iElement = 0; iElement < sub->vMaps.size(); iElement++) {
      string elementName = sub->vMaps[iElement]->getName();
      ConstMapIter j = mapElements.find(elementName);
      if (j==mapElements.end()) 
	throw PhotometryError("PhotoMapCollection::rebuildParameterVector could not find"
			      " map element with name " + elementName);
      if (!j->second.atom)
	throw PhotometryError("PhotoMapCollection element " + elementName + " is not atomic"
			      " in rebuildParameterVector.");
      sub->vStartIndices[iElement] = j->second.startIndex;
      sub->vNSubParams[iElement] = j->second.nParams;
      sub->vMapNumbers[iElement] = j->second.number;
    }
    sub->countFreeParameters();
  }
}


//////////////////////////////////////////////////////////////
// Member maintenance
//////////////////////////////////////////////////////////////

vector<string>
PhotoMapCollection::allMapNames() const {
  vector<string> output;
  for (ConstMapIter i = mapElements.begin();
       i != mapElements.end();
       ++i) 
    output.push_back(i->first);
  return output;
}

void 
PhotoMapCollection::learnMap(const PhotoMap& pm,
			     bool duplicateNamesAreExceptions,
			     bool rebuildIndices) {
  if (mapExists(pm.getName())) {
    if (duplicateNamesAreExceptions)
      throw PhotometryError("Duplicate map name in PhotoMapCollection::learnMap at "
			    + pm.getName());
    // If not throwing an exception, we will just ignore this duplicate-named map.
    return;
  }
  const SubMap* sm = dynamic_cast<const SubMap*> (&pm);
  if (sm) {
    if (sm->vMaps.size()==1 && pm.getName()==sm->vMaps.front()->getName()) {
      // If this is a single-element SubMap that has same name as its dependent, then we'll
      // just learn the dependent:
      learnMap(*sm->vMaps.front(), duplicateNamesAreExceptions, rebuildIndices);
      return;
    } 
    // create a new compound map from this, learning each dependency:
    MapElement mel;
    for (int i = 0; i< sm->vMaps.size(); i++) {
      if (sm->getName() == sm->vMaps[i]->getName())
	throw PhotometryError("PhotoMapCollection::learnMap encountered SubMap that has "
			      "a dependence on PhotoMap with the same name: " + sm->getName());
      mel.subordinateMaps.push_back(sm->vMaps[i]->getName());
      learnMap(*sm->vMaps[i], duplicateNamesAreExceptions, rebuildIndices);
    }
    // Register this compound map
    mapElements.insert(std::pair<string,MapElement>(pm.getName(), mel));

  } else {
    //  This should be an atomic PhotoMap; simply add it to our element list:
    MapElement mel;
    mel.atom = pm.duplicate();
    mel.isFixed = mel.atom->nParams()==0;
    mapElements.insert(std::pair<string,MapElement>(pm.getName(), mel));
  }

  // Check that new map did not introduce circularity (though I do not know how it could)
  checkCircularDependence(pm.getName());

  // re-index everything
  if (rebuildIndices) rebuildParameterVector();
}

void
PhotoMapCollection::learn(PhotoMapCollection& rhs, bool duplicateNamesAreExceptions) {
  if (&rhs == this) return;	// No need to learn self
  for (ConstMapIter iMap = rhs.mapElements.begin();
       iMap != rhs.mapElements.end();
       ++iMap) {
    const MapElement& incoming = iMap->second;
    if (mapExists(iMap->first)) {
      // incoming map duplicates and existing name.
      if (duplicateNamesAreExceptions) 
	throw PhotometryError("learn(PhotoMapCollection) with duplicate map name "
			      + iMap->first);
      // Duplicate will just be ignored. 
    } else {
      // A new mapName for us.  Add its mapElement to our list.
      MapElement mel;
      if (iMap->second.atom) mel.atom = iMap->second.atom->duplicate();
      mel.isFixed = iMap->second.isFixed;
      mel.subordinateMaps = iMap->second.subordinateMaps;
      mapElements.insert(std::pair<string,MapElement>(iMap->first, mel));
      // Note cannot introduce circularity if the top-level map is new.
    }
  } // end input map loop

  // And re-index everything
  rebuildParameterVector();
}

// Define a new PhotoMap that is compounding of a list of other PhotoMaps.  Order
// of the list is order of transform from mapIn to mapOut
void 
PhotoMapCollection::defineChain(string chainName, const list<string>& elements) {
  if (mapExists(chainName)) 
    throw PhotometryError("PhotoMapCollection::defineChain with duplicate name: " 
			  + chainName);
  // Check that elements exist
  for (list<string>::const_iterator i = elements.begin();
       i != elements.end();
       ++i) 
    if (!mapExists(*i))
      throw PhotometryError("PhotoMapCorrection::defineChain with unknown photo map element: "
			    + *i);
  mapElements.insert(std::pair<string, MapElement>(chainName, MapElement()));
  mapElements[chainName].subordinateMaps = elements;
  // Note that adding a new chain does not change parameter vector assignments
  // Nor can it introduce circular dependence if chain name is new and all elements exist.
}

// Return pointer to a SubMap realizing the named magnitude transformation
SubMap* 
PhotoMapCollection::issueMap(string mapName) {
  if (!mapExists(mapName))
    throw PhotometryError("PhotoMapCollection::issueMap requested for unknown PhotoMap: "
			  +  mapName);
  MapElement& el = mapElements[mapName];

  if (!el.realization) {
    // Create a realization if one does not exist
    list<string> atomList = orderAtoms(mapName);
    vector<int> startIndices(atomList.size(), 0);
    vector<int> nParams(atomList.size(), 0);
    vector<int> mapNumbers(atomList.size(), -1);

    list<PhotoMap*> atoms;
    int index=0;
    for (list<string>::const_iterator i = atomList.begin();
	 i != atomList.end();
	 ++i, ++index) {
      Assert(mapElements[*i].atom);	// All elements should be atomic
      atoms.push_back(mapElements[*i].atom);
      // fill in its indices into master vector:
      startIndices[index] = mapElements[*i].startIndex;
      nParams[index] = mapElements[*i].isFixed ? 0 : mapElements[*i].nParams;
      mapNumbers[index] = mapElements[*i].number;
    }
    SubMap* sm = new SubMap(atoms, mapName, true);
    sm->vStartIndices = startIndices;
    sm->vNSubParams = nParams;
    sm->vMapNumbers = mapNumbers;
    sm->countFreeParameters();
    el.realization = sm;
  }
  return el.realization;
}

PhotoMap*
PhotoMapCollection::cloneMap(string mapName) const {
  if (!mapExists(mapName))
    throw PhotometryError("PhotoMapCollection::issueMap requested for unknown PhotoMap: "
			  +  mapName);
  const MapElement& el = mapElements.find(mapName)->second;

  if (el.atom) 
    return el.atom->duplicate();

  // A composite one we will create as a SubMap:
  list<PhotoMap*> vMaps;
  for (list<string>::const_iterator i = el.subordinateMaps.begin();
	 i != el.subordinateMaps.end();
	 ++i) 
    vMaps.push_back(cloneMap(*i));
  SubMap* retval = new SubMap(vMaps, mapName, false);
  // Clean up the stray clones since they've been duplicated by SubMap:
  for (list<PhotoMap*>::iterator i = vMaps.begin();
       i != vMaps.end();
       ++i)
    delete *i;
  return retval;
}

set<string>
PhotoMapCollection::dependencies(string target) const {
  // Find target MapElement and add target to output dependency list.  
  // Assume no circular dependence.
  ConstMapIter iTarget = mapElements.find(target);
  if (iTarget == mapElements.end()) 
    throw PhotometryError("PhotoMapCollection has no element " + target + " in dependencies()");
  set<string> output;
  output.insert(target);

  // Call this routine recursively for all the map elements that the target depends on
  for (list<string>::const_iterator i = iTarget->second.subordinateMaps.begin();
       i != iTarget->second.subordinateMaps.end();
       ++i) {
    set<string> subs = dependencies(*i);
    output.insert(subs.begin(), subs.end());
  }

  return output;
}

bool
PhotoMapCollection::dependsOn(const string mapName, const string targetName) const {
  auto melptr = mapElements.find(mapName);
  if (melptr == mapElements.end()) 
    throw PhotometryError("PhotoMapCollection has no element " + mapName + " in dependsOn()");

  if (mapName==targetName) return true;  // True if map is the target

  // Call this routine recursively for all the map elements that the target depends on
  for (auto& i : melptr->second.subordinateMaps)
    if (dependsOn(i, targetName)) return true;

  return false;
}


// Return a list of the atomic elements of the named map, in order of
// their application to data.  Assumes no circular dependence.
list<string> 
PhotoMapCollection::orderAtoms(string mapName) const {
  if (!mapExists(mapName))
    throw PhotometryError("PhotoMapCollection::orderAtoms requested for unknown map: "
			  + mapName);
  list<string> retval;
  const MapElement& m = mapElements.find(mapName)->second;
  if (m.atom) {
    // Simple atomic element just returns itself
    retval.push_back(mapName);
    return retval;
  }

  // Otherwise we have a compound map at work here.
  // Call recursively to all subordinateMaps:
  for (list<string>::const_iterator i = m.subordinateMaps.begin();
       i != m.subordinateMaps.end();
       ++i) {
    list<string> subList = orderAtoms(*i);
    retval.splice(retval.end(), subList);
  }
  return retval;
}

void
PhotoMapCollection::checkCircularDependence(string mapName,
					    const set<string>& ancestors) const {
  if (!mapExists(mapName))
    throw PhotometryError("Unknown mapName in PhotoMapCollection::checkCircularDependency: "
			  + mapName);
  if (ancestors.count(mapName))
    throw PhotometryError("Circular dependency in PhotoMapCollection at map "
			  + mapName);
  // Then call recursively to all subordinateMaps, adding self to ancestor list:
  set<string> ancestorsPlusSelf(ancestors);
  ancestorsPlusSelf.insert(mapName);
  const MapElement& m = mapElements.find(mapName)->second;
  for (list<string>::const_iterator i = m.subordinateMaps.begin();
       i != m.subordinateMaps.end();
       ++i) 
    checkCircularDependence(*i, ancestorsPlusSelf);
}

void
PhotoMapCollection::checkCompleteness() const {
  // Loop through all the (non-atomic) maps, throw exception if the refer to 
  // non-existent map elements or to themselves
  for (ConstMapIter i = mapElements.begin();
       i != mapElements.end();
       ++i)
    checkCircularDependence(i->first);
}


//////////////////////////////////////////////////////////////
// YAML (De-) Serialization
//////////////////////////////////////////////////////////////

void
PhotoMapCollection::writeSingleMap(YAML::Emitter& os, const MapElement& mel, string name)  const {
  if (mel.atom) {
    // Atomic map, give all details:
    const ColorTerm* ct = dynamic_cast<const ColorTerm*> (mel.atom);
    if (ct) {
      // We have a color term which is wrapping another map.
      // Make a wrapper node of Type=Color and write the wrapped function
      os << YAML::Key << name
	 << YAML::Value
	 << YAML::BeginMap
	 << YAML::Key << "Type" << YAML::Value << ColorTerm::type()
	 << YAML::Key << "Reference" << YAML::Value << ct->refColor()
	 << YAML::Key << "Function" << YAML::Value;
      ct->map()->write(os);
      os << YAML::EndMap;
    } else {
      // A simple atomic node.  Let it write itself
      os << YAML::Key << name
	 << YAML::Value;
      mel.atom->write(os);
    }
  } else {
    // For compound map, just give names of constituents
    os << YAML::Key << name
       << YAML::Value 
       << YAML::BeginMap
       << YAML::Key << "Type"
       << YAML::Value << SubMap::type()
       << YAML::Key << "Elements"
       << YAML::BeginSeq;
    for (auto subname : mel.subordinateMaps) {
      os << subname;
    }
    os << YAML::EndSeq
       << YAML::EndMap;
  }
}

void 
PhotoMapCollection::write(ostream& os, string comment) const {
  // Create a YAML emitter.  Produce a map where one
  // key is the magic word, and the rest are names of maps.
  YAML::Emitter out;
  out << YAML::BeginMap;
  if (comment.size()>0)
    out << YAML::Comment(comment);
  out << YAML::Key << magicKey
      << YAML::Value << "This is a serialized PhotoMapCollection";
  for (auto& melpair : mapElements) {
    writeSingleMap(out, melpair.second, melpair.first);
  }
  out << YAML::EndMap;
  // Send the YAML to the output stream
  os << out.c_str() << endl;
}

void 
PhotoMapCollection::writeMap(ostream& os, string name, string comment) const {
  if (!mapExists(name))
    throw PhotometryError("PhotoMapCollection::writeMap() for unknown map: " + name);
  set<string> allmaps = dependencies(name);
  YAML::Emitter out;
  out << YAML::BeginMap;
  if (comment.size()>0)
    out << YAML::Comment(comment);
  out << YAML::Key << magicKey
      << YAML::Value << "This is a serialized PhotoMapCollection";
  // Write the map elements
  for (auto mapname : allmaps) {
    auto j = mapElements.find(mapname);
    writeSingleMap(out, j->second, mapname);
  }
  out << YAML::EndMap;
  // Send the YAML to the output stream
  os << out.c_str() << endl;
}

bool
PhotoMapCollection::read(istream& is, string namePrefix) {
  string buffer;
  try {
    YAML::Node root = YAML::Load(is);
    if (!root.IsMap() || !root[magicKey]) {
      // A valid YAML document that does not have the structure 
      // and key that our files should have.
      return false;
    }

    // Every entry in the root node map is a map to enter
    for (YAML::const_iterator i = root.begin();
	 i != root.end();
	 ++i) {
      // Skip the map key that is our magic word
      string keyName = i->first.as<string>();
      if (keyName == magicKey) continue;

      string name = namePrefix + keyName;
      const YAML::Node& node=i->second;

      if (!node.IsMap() || !node["Type"]) 
	throw PhotometryError("Non-map or missing Type key in PhotoMapCollection YAML entry " +
			      name);
      string mapType = node["Type"].as<string>();
      if (stringstuff::nocaseEqual(mapType,ColorTerm::type())) {
	// A color term has the function it modifies in a child node:
	if (!node["Function"])
	throw PhotometryError("Missing Function key in ColorTerm YAML entry " +
			      name);
	double reference = 0.; // Default reference color
	if (node["Reference"])
	  reference = node["Reference"].as<double>();
	PhotoMap* pm = createAtomFromNode(node["Function"], "");

	// create the MapElement
	MapElement mel;
	mel.atom = new ColorTerm(pm, reference, name);
	mapElements.insert(std::pair<string,MapElement> (name, mel));
	
      } else if (stringstuff::nocaseEqual(mapType,SubMap::type())) {
	// Composite maps: all names stored in a child node
	if (!node["Elements"])
	  throw PhotometryError("Missing Elements key in SubMap YAML entry " + name);
	list<string> submaps = node["Elements"].as<list<string> >();

	// Create the MapElement
	MapElement mel;
	mel.subordinateMaps = submaps;
	mapElements.insert(std::pair<string,MapElement>(name, mel));

      } else {
	// Should have some kind of atomic map.  Call a function to parse the node.
	// This will throw exception if we have an unknown Type of map
	
	// create the MapElement
	MapElement mel;
	mel.atom = createAtomFromNode(node,name);
	mel.isFixed = mel.atom->nParams()==0;  // Set fixed if no free params
	mapElements.insert(std::pair<string,MapElement> (name, mel));
      }

    } // end loop through root node map members

  } catch (YAML::Exception& e) {
    /**/cerr << "PhotoMapCollection::read() caught: " << e.what() << endl;
    // Return false if not a valid YAML file
    return false;
  }

  // Check all maps for circular dependence and completeness.
  checkCompleteness();

  // Recalculate indices
  rebuildParameterVector();

  return true;
}

PhotoMap*
PhotoMapCollection::createAtomFromNode(const YAML::Node& node, string name) const {
  // Access the static arrays of map type names and creators
  if (!node["Type"])
    throw PhotometryError("Missing Type key in createAtomFromNode at " + name);

  string mapType = node["Type"].as<string>();
  auto it = mapCreators.find(mapType);
  if (it==mapCreators.end())
    throw PhotometryError("PhotoMapCollection does not recognize PhotoMap Type " + mapType);
  return (*(it->second))(node, name);
}

//////////////////////////////////////////////////////////////////////////////

string
PhotoMapCollection::atomHavingParameter(int parameterIndex) const {
  // Find the atomic map that this parameter number influences
  for (auto& melpair : mapElements) {
    if (!melpair.second.atom) continue;
    if ( parameterIndex >= melpair.second.startIndex &&
	 parameterIndex < melpair.second.startIndex + melpair.second.nParams)
      return melpair.first;
  }
  return "";	// Nothing in range, return empty string.
}

void
PhotoMapCollection::parameterIndicesOf(string mapname,
				       int& startIndex,
				       int& nParams) const {
  if (!mapExists(mapname))
    throw PhotometryError("PhotoMapCollection::parameterIndicesof() unknown map: "
			  + mapname);
  startIndex = 0;
  nParams = 0;
  auto& mel = mapElements.at(mapname);
  if (!mel.atom || mel.isFixed) return;
  startIndex = mel.startIndex;
  nParams = mel.nParams;
  return;
}
