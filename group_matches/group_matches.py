# Function usage tree
#	
# 	/--	get_doc_condensed_matches
# 		|
# 		+--	get_matches
# 		|	|
# 		|	+--	rec_fetch_parent
# 		|
# 		+--	group_by_overlap_and_ontology
# 		|	|
# 		|	+--	exact_path_grouping
# 		|	|
# 		|	+--	detect_path_overlaps
# 		|	|
# 		|	+--	condense_group
# 		|		|
# 		|		+--	group_codes
# 		|		|	|
# 		|		|	+--	get_equivalent_links
# 		|		|	|
# 		|		|	+--	get_subclass_links
# 		|		|	|
# 		|		|	+--	group_by_links
# 		|		|
# 		|		+--	pick_top_matches
# 		|		|	|
# 		|		|	+--	group_by_links
# 		|		|	|
# 		|		|	+--	assign_group_distances
# 		|		|		|
# 		|		|		+--	rec_assign_group_distances
# 		|		|
# 		|		+--	add_nearest_efo
# 		|		|	|
# 		|		|	+--	find_closest_efo
# 		|		|		|
# 		|		|		+--	get_eq
# 		|		|
# 		|		+--	label_groups
# 		|		|	|
# 		|		|	+--	get_top_label
# 		|		|		|
# 		|		|		+--	get_labels
# 		|		|
# 		|		+--	efo_grouping
# 		|			|
# 		|			+--	group_by_links
# 		|
# 		+--	matches_grouping
# 			|
# 			+--	group_matches_by_links

import itertools as it

# join the dicts preserving the lower distance (v)
def join_entities_dict(source_dict, dest_dict):
    for k,v in source_dict.items():
        if k in dest_dict:
            if dest_dict[k] > v: dest_dict[k] = v 
        else: dest_dict[k] = v
    
    return dest_dict

def rec_fetch_parent(word_id, sentence):
    parent_ids = set()
    word = sentence['words'][word_id]
    if 'parent' in word.keys():
        for parent_id in word['parent']:
            if parent_id == word_id: parent_ids.add(word_id) # break the infinite loop
            else: parent_ids.update(rec_fetch_parent(parent_id, sentence))
    else: parent_ids.add(word_id)
    
    return parent_ids

def get_equivalent_links(onto_codes, equivalence_index):
    links = set()
    for code1,code2 in it.combinations(onto_codes,2):
        try: groups_1 = equivalence_index[code1[0]] | {code1[0]}
        except: continue

        try: groups_2 = equivalence_index[code2[0]] | {code2[0]}
        except: continue

        if groups_1 & groups_2:
            links.add(tuple(sorted([code1, code2])))
    
    return links

def get_subclass_links(onto_codes, equivalence_index, equivalence_index_r, subclass_index):
    links = set()
    for code1,code2 in it.combinations(onto_codes,2):
        try: code1_entities = equivalence_index_r[str(equivalence_index[code1[0]])] | {code1[0]}
        except: code1_entities = {code1[0]}
        try: code2_entities = equivalence_index_r[str(equivalence_index[code2[0]])] | {code2[0]}
        except: code2_entities = {code2[0]}
        
        code1_subclass_entities = code1_entities.copy()
        for e in code1_entities:
            e_source, e_code = e.split(':')
            try: code1_subclass_entities.update(subclass_index[e_source][e])
            except: continue
        
        code2_subclass_entities = code2_entities.copy()
        for e in code2_entities:
            e_source, e_code = e.split(':')
            try: code2_subclass_entities.update(subclass_index[e_source][e])
            except: continue
                
        if (code1_subclass_entities & code2_entities) | (code2_subclass_entities & code1_entities):
            links.add((code1, code2)) 
    
    return links

def group_by_links(links, groups):
    links = [set(l) for l in links]
    
    for link in links:
        join_groups = set()
        for i,g in enumerate(groups):
            if g & link:
                join_groups.add(i)

        if join_groups:
            merge_group_i = sorted(join_groups, reverse=True)[0]

            for i in sorted(join_groups, reverse=True):
                if i == merge_group_i: continue
                groups[merge_group_i].update(groups[i])
                groups[i] = set()

    groups = [g for g in groups if len(g)]

    return groups

def group_codes(onto_codes, equivalence_index, equivalence_index_r, subclass_index):
    groups = [{code} for code in onto_codes]
    equivalent_links = get_equivalent_links(onto_codes, equivalence_index)
    subclass_links = get_subclass_links(onto_codes, equivalence_index, equivalence_index_r, subclass_index)
    
    groups_dict = {i:v for i,v in enumerate(group_by_links(equivalent_links, groups))}

    groups_dict_r = {}
    for k,vs in groups_dict.items():
        for v in vs:
            groups_dict_r[v] = k

    group_rels = set()
    for c_start, c_end in subclass_links:
        c_start = groups_dict_r[c_start]
        c_end = groups_dict_r[c_end]

        if c_start == c_end: continue

        group_rels.add((c_start, '-[subclass_of]->', c_end))

    group_rels
    
    return groups_dict, group_rels

def rec_assign_group_distances(group, links_dict, distance, group_distances):
    if group in group_distances.keys():
        if group_distances[group] > distance:
            group_distances[group] = distance
    else:
        group_distances[group] = distance
        
    if group in links_dict.keys():
        for g in links_dict[group]:
            distances = rec_assign_group_distances(g, links_dict, distance+1, group_distances)
            for k,v in distances.items():
                if k in group_distances.keys():
                    if group_distances[k] > v:
                        group_distances[k] = v
                else:
                    group_distances[k] = v
                
    return group_distances

def assign_group_distances(links):
    links_dict = {}
    for l in links:
        try: links_dict[l[0]].add(l[1])
        except: links_dict[l[0]] = {l[1]}

    links_dict_r = {}
    for k,vs in links_dict.items():
        for v in vs:
            try: links_dict_r[v].add(k)
            except: links_dict_r[v] = {k}

    roots = set(links_dict.keys()) - set(links_dict_r.keys())

    group_distances = {}
    for root in roots:
        distances = rec_assign_group_distances(root, links_dict, 0, {})
        for k,v in distances.items():
            if k in group_distances.keys():
                if group_distances[k] > v:
                    group_distances[k] = v
            else:
                group_distances[k] = v

    return group_distances

def pick_top_matches(groups_dict, group_rels, source_ranks, predicate_type_ranks):
    group_groups = group_by_links({(r[0],r[2]) for r in group_rels}, [{g} for g in groups_dict.keys()])
    
    group_picked_matches = []
    for group_group in group_groups:
        # dynamically pick the groups here (not rely on precomputed lookup tables)
        group_group_rels = [(g[0],g[2]) for g in group_rels if ({g[0],g[2]} & group_group)]

        efo_groups = set()
        for group in group_group:
            g_sources = {m[0].split(':')[0] for m in groups_dict[group]}
            if 'efo' in g_sources:
                efo_groups.add(group)

        if len(group_group_rels):
            group_distances = assign_group_distances(group_group_rels)
            try: max_distance = max(group_distances.values())
            except: max_distance = 0

            top_groups = set()
            for k,v in group_distances.items():
                if v >= max_distance:
                    top_groups.add(k)

            try: efo_groups = {sorted(efo_groups, key=lambda x:group_distances[x], reverse=True)[0]}
            except: efo_groups = set()
        else:
            top_groups = group_group.copy()


        picked_matches = set()
        for group in top_groups:
            group_matches = groups_dict[group]

            top_predicate_rank = max([predicate_type_ranks[m[1]] for m in group_matches])
            top_predicate_matches = [m for m in group_matches if predicate_type_ranks[m[1]] == top_predicate_rank]

            top_source_rank = max([source_ranks[m[0].split(':')[0]] for m in top_predicate_matches])
            top_matches = [m for m in top_predicate_matches if source_ranks[m[0].split(':')[0]] == top_source_rank]

            picked_matches.update({m for m in top_matches})

        for group in efo_groups:
            group_matches = groups_dict[group]

            efo_matches = {m for m in group_matches if m[0].split(':')[0] == 'efo'}
            picked_matches.update(efo_matches)

        group_picked_matches.append(picked_matches)

    return group_picked_matches

def get_eq(code, eq_groups_index, eq_groups_index_r):  # Get equivalence group ID for the code, then return the codes in that group
    eq_codes = {code}
    try: 
        eq_groups = eq_groups_index_r[code]
        for eq_group in eq_groups:
            eq_codes.update(eq_groups_index[eq_group])
    except: pass
    
    return eq_codes

# returns all related EFO IDs with their distances, negative distances and superclasses, positive distances are subclasses (superclasses are better)
def find_closest_efo(subject_code, eq_groups_index, eq_groups_index_r, subclass_indexes, rev_subclass_indexes):
    all_related_codes = set()
    
    eq_codes = get_eq(subject_code, eq_groups_index, eq_groups_index_r)  # get equivalent terms of subject_code
    
    # get all subclass and superclasses of the equivalent codes
    for code in eq_codes:
        all_related_codes.add((code, code.split(':')[0], 0))
        for source,subclass_index in subclass_indexes.items():
            try: all_related_codes.update({(v, v.split(':')[0], d) for v,d in subclass_index[code].items()})  # subclasses
            except: pass
            try: all_related_codes.update({(v, v.split(':')[0], -d) for v,d in rev_subclass_index[code].items()})  # superclasses
            except: pass
    
    # get equivalent codes of all of the subclasses and superclasses
    all_related_codes_expanded = set()
    for code,source,distance in all_related_codes:
        all_related_codes_expanded.add((code,source,distance))
        all_related_codes_expanded.update((c, c.split(':')[0], distance) for c in get_eq(code, eq_groups_index, eq_groups_index_r))
    
    # filter just EFO IDs and sort by distance (shortest distance first)
    sorted_efo_related_codes = sorted([(c,d) for c,s,d in all_related_codes_expanded if s=='efo'], key=lambda x:abs(x[1]))
    return sorted_efo_related_codes

def get_labels(index_conn, onto_id):
    results = index_conn.execute(f"select string, predicate_type from strings where onto_id='{onto_id}'")
    labels = {}
    for string, predicate_type in results:
        try: labels[predicate_type].add(string)
        except: labels[predicate_type] = {string}
    
    return labels

def get_top_label(index_conn, onto_id, predicate_type_ranks):
    labels = get_labels(index_conn, onto_id)
    if len(labels):
        top_predicate_rank = max([predicate_type_ranks[k] for k in labels.keys()])
        top_predicates = {l for k,v in labels.items() for l in v if predicate_type_ranks[k] == top_predicate_rank}
        return sorted(top_predicates, key=lambda x:len(x))[0]
    else: return None
    
def label_groups(picked_groups, index_conn, predicate_type_ranks):
    picked_labelled_groups = []
    for g in picked_groups:
        picked_labelled_groups.append({(m[0], m[1], get_top_label(index_conn, m[0], predicate_type_ranks)) for m in g})
        
    return picked_labelled_groups

def add_nearest_efo(groups, equivalent_entities_groups_index, equivalent_entities_groups_index_r, disease_hierarchy_distance_index, rev_disease_hierarchy_distance_index):
    for i,group in enumerate(groups):
        efo_ids = set()
        for c,p in group:
            if c.split(':')[0] == 'efo':
                efo_ids.add((c,0))
        if len(efo_ids) == 0:
            for c,p in group:
                code_efo_ids = find_closest_efo(c, equivalent_entities_groups_index, equivalent_entities_groups_index_r, disease_hierarchy_distance_index, rev_disease_hierarchy_distance_index)
                try: 
                    min_d = min([abs(d) for c,d in code_efo_ids if d<=0])
                    efo_ids.update({(c,d) for c,d in code_efo_ids if abs(d)==min_d})
                except: pass
        
        try: 
            min_d = min([abs(d) for c,d in efo_ids])
            groups[i].update({(c,f'nearest_efo:{d}') for c,d in efo_ids if abs(d)==min_d})
        except: pass
    
    return groups

def efo_grouping(groups):  # group based on common EFO IDs if the distance is 0
    links = set()
    for g1,g2 in it.combinations(range(len(groups)),2):
        group1 = {(c,l) for c,p,l in groups[g1] if (not p.split(':')[0]=='nearest_efo') or (p=='nearest_efo:0')}
        group2 = {(c,l) for c,p,l in groups[g2] if (not p.split(':')[0]=='nearest_efo') or (p=='nearest_efo:0')}
        if group1 & group2:
            links.add(sorted([g1,g2]))

    new_groups = group_by_links(links, groups)

    return new_groups

def condense_group(group, index_conn, equivalent_entities_groups_index, equivalent_entities_groups_index_r, disease_hierarchy_index, disease_hierarchy_distance_index, rev_disease_hierarchy_distance_index, source_ranks, predicate_type_ranks):
    groups_dict, group_rels = group_codes(group, equivalent_entities_groups_index_r, equivalent_entities_groups_index, disease_hierarchy_index)
    picked_groups = pick_top_matches(groups_dict, group_rels, source_ranks, predicate_type_ranks)
    picked_groups_w_efo = add_nearest_efo(picked_groups, equivalent_entities_groups_index, equivalent_entities_groups_index_r, disease_hierarchy_distance_index, rev_disease_hierarchy_distance_index)
    picked_labelled_groups_w_efo = label_groups(picked_groups_w_efo, index_conn, predicate_type_ranks)
    condensed_groups = efo_grouping(picked_labelled_groups_w_efo)
    
#     for i,group in enumerate(condensed_groups):
#         condensed_groups[i] = {(c,l) for c,p,l in group}
    
    return condensed_groups

# def get_matches(sentence_id, indi_conn):
#     results = indi_conn.execute(f'select id, onto_id, predicate_type, source, locs from matches where sentence_id={sentence_id}')
#     expanded_sentence = eval(list(indi_conn.execute(f"select expanded_sentence from sentences where id={sentence_id}"))[0][0])
    
#     matches = set()
#     for match_id, onto_id, predicate_type, source, locs in results:
#         locs = eval(locs)
#         original_locs = set()
#         for path in locs:
#             original_path = set()
#             for l in path:
#                 original_path.update(rec_fetch_parent(l, expanded_sentence))
#             original_locs.add(tuple(original_path))
#         matches.add((onto_id, predicate_type, tuple(original_locs)))
    
#     return matches

def exact_path_grouping(matches):
    exact_path_groups = {}
    for onto_id, predicate_type, paths in matches:
        for p in paths:
            try: exact_path_groups[p].add((onto_id,predicate_type))
            except: exact_path_groups[p]= {(onto_id,predicate_type)}
    
    return exact_path_groups

def detect_path_overlaps(exact_path_groups):
    overlaps = set()
    for p1,p2 in it.combinations(exact_path_groups.keys(),2):
        p1_set = set(p1)
        p2_set = set(p2)
        if len(p1_set - p2_set) == 0:  # p1 subset of p2
            overlaps.add((p1,'subset',p2))
            continue
        if len(p2_set - p1_set) == 0:  # p2 subset of p1
            overlaps.add((p2,'subset',p1))
            continue
        if p2_set & p1_set:  # overlap
            overlaps.add((p1,'overlap',p2))
            continue
            
    return overlaps

def group_by_overlap_and_ontology(matches, index_conn, equivalent_entities_groups_index, equivalent_entities_groups_index_r, disease_hierarchy_index, disease_hierarchy_distance_index, rev_disease_hierarchy_distance_index, source_ranks, predicate_type_ranks):
    exact_path_groups = exact_path_grouping(matches)
    match_path_overlaps = detect_path_overlaps(exact_path_groups)
    
    superset_paths = {p[0] for p in match_path_overlaps if p[1] == 'subset'}
    
    condensed_matches = {}
    for path, group in exact_path_groups.items():
        if path in superset_paths: continue
        condensed_matches[path] = condense_group(group, index_conn, equivalent_entities_groups_index, equivalent_entities_groups_index_r, disease_hierarchy_index, disease_hierarchy_distance_index, rev_disease_hierarchy_distance_index, source_ranks, predicate_type_ranks)
        
    return condensed_matches

def matches_grouping(groups, source_ranks):  # group based on common EFO IDs if the distance is 0
    groups = list(groups.items())
    links = set()
    for g1,g2 in it.combinations(range(len(groups)),2):
        group1 = {c for c,p,l in groups[g1][0] if (not p.split(':')[0]=='nearest_efo') or (p=='nearest_efo:0')}
        group2 = {c for c,p,l in groups[g2][0] if (not p.split(':')[0]=='nearest_efo') or (p=='nearest_efo:0')}
        if group1 & group2:
            links.add(tuple(sorted([g1,g2])))
    
    new_groups = group_matches_by_links(links, [{(tuple(k),tuple(v))} for k,v in groups])
    
    new_groups_dict = {}
    for s in new_groups:
        k = set()
        v = set()
        for s_k,s_v in s:
            k.update(set(s_k))
            v.update(set(s_v))

        matches = {}
        for c,p,l in k:
            p = 'direct' if ((not p.split(':')[0]=='nearest_efo') or (p=='nearest_efo:0')) else 'distant'
            if (c,l) in matches:
                if p == 'direct':
                    matches[(c,l)] = p
            else: matches[(c,l)] = p

        matches = {(c,p,l) for (c,p),l in matches.items()}
        matches = tuple(sorted(matches, key=lambda x:(source_ranks[x[0].split(':')[0]], 0 if x[2] is None else 1, x[0]), reverse=True))
        new_groups_dict[matches] = v
    
    return new_groups_dict

def group_matches_by_links(links, groups):
    links = [set(l) for l in links]
    
    for link in links:
        join_groups = set()
        for i,g in enumerate(groups):
            if {i} & link:
                join_groups.add(i)

        if join_groups:
            merge_group_i = sorted(join_groups, reverse=True)[0]

            for i in sorted(join_groups, reverse=True):
                if i == merge_group_i: continue
                groups[merge_group_i].update(groups[i])
                groups[i] = set()

    groups = [g for g in groups if len(g)]

    return groups

def sort_score_f(x):
    try: 
        predicate_score = predicate_type_ranks[x[1]]
    except:
        if x[1] == "nearest_efo:0":
            predicate_score = 7
        else:
            predicate_score = 0

    source_score = source_ranks[x[0].split(':')[0]]

    name_score = 0 if x[2] is None else 1

    return (predicate_score, source_score, name_score, x[0])

def get_doc_condensed_matches(  sentences, 
                                indi_conn, 
                                index_conn, 
                                equivalent_entities_groups_index, 
                                equivalent_entities_groups_index_r, 
                                disease_hierarchy_index, 
                                disease_hierarchy_distance_index, 
                                rev_disease_hierarchy_distance_index, 
                                source_ranks =         {'icd11': 0,
                                                        'doid': 1,
                                                        'orphanet': 1,
                                                        'mondo': 1,
                                                        'icd10': 2,
                                                        'snomed': 3,
                                                        'mesh': 4,
                                                        'efo': 5}, 
                                predicate_type_ranks = {'note': 0,
                                                        'desc': 1,
                                                        'oboInOwl:hasBroadSynonym': 2,
                                                        'oboInOwl:hasNarrowSynonym': 2,
                                                        'oboInOwl:hasRelatedSynonym': 3,
                                                        'umls': 3,
                                                        'synonym': 4,
                                                        'efo:alternative_term': 4,
                                                        'altLabel': 4,
                                                        'oboInOwl:hasExactSynonym': 5,
                                                        'label': 6,
                                                        'rdfs:label': 6,
                                                        'prefLabel': 7,
                                                        'skos:prefLabel': 7,
                                                        'pref': 7,
                                                        'title': 7}):
    
    doc_condensed_matches = {}
    for s_id, loc, sentence, expanded_sentence, matches in sentences:
        # matches = get_matches(s_id, indi_conn)
        condensed_matches = group_by_overlap_and_ontology(matches, index_conn, equivalent_entities_groups_index, equivalent_entities_groups_index_r, disease_hierarchy_index, disease_hierarchy_distance_index, rev_disease_hierarchy_distance_index, source_ranks, predicate_type_ranks)
        for path,group in condensed_matches.items():
            doc_condensed_matches[(loc,path)] = [tuple(sorted(m, key=lambda x:sort_score_f(x), reverse=True)) for m in group]

    doc_condensed_matches_r = {}
    for k,vs in doc_condensed_matches.items():
        matches = []
        for v in vs:
            try: doc_condensed_matches_r[matches].add(k)
            except: doc_condensed_matches_r[matches] = {k}
    
    doc_condensed_matches_r = matches_grouping(doc_condensed_matches_r, source_ranks)
    
    doc_condensed_matches_reformated = []
    for k,v in doc_condensed_matches.items():
        doc_condensed_matches_reformated.append({'loc': {'sentence_loc': k[0], 'path': k[1]}, 'match_groups': {i: list(matches) for i,matches in enumerate(v)}})
    
    doc_condensed_matches_r_reformated = []
    for k,vs in doc_condensed_matches_r.items():
        doc_condensed_matches_r_reformated.append({'locs': [{'sentence_loc': v[0], 'path': v[1]} for v in vs], 'matches': list(k)})
    
    return doc_condensed_matches_reformated, doc_condensed_matches_r_reformated
