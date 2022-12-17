import _ from "lodash";

import modelFeatures from "../feature_data";
import { recursiveMerge } from "./utils";

/**
 * Filter nodes based on path and modify node data.
 * @param {Array} nodes - Array of nodes to be filtered.
 * @param {string[]} paths - Array of node paths
 * @param {Object[]} pathData - Array of objects containing path and node data.
 * @returns {Object[]} - The filtered tree (with possibly modified data).
 */
export function filterTree(nodes, paths, pathData) {
    return nodes.reduce((acc, node) => {
        if (paths.includes(node.label)) {
            let modified = {};
            let data = pathData && pathData.find((item) => item.path === node.label)?.data;
            if (data) {
                data = recursiveMerge(node.data, data);
                modified = { ...node, data };
            }
            if (node.children?.length) {
                const children = filterTree(node.children, paths, pathData);
                if (children.length) modified = { ...node, ...modified, children };
            }
            // eslint-disable-next-line no-unused-expressions
            _.isEmpty(modified) ? acc.push(node) : acc.push(modified);
        }
        return acc;
    }, []);
}

/**
 * Selects a subset of the model tree for a given application
 * @param {Object|Array} tree - Tree object or array of nodes
 * @param {string} appName
 * @param {string} version
 * @param {string} executable
 * @param {string} build
 */
export function getModelTreeByApplication({
    tree,
    appName,
    version,
    executable,
    build = "default",
}) {
    const featureSubtree = modelFeatures[appName][version][executable][build];
    const nodePaths = featureSubtree.map((item) => item.path);
    const pathData = featureSubtree.filter((item) => "data" in item);
    const nodes = Array.isArray(tree) ? tree : [tree];
    return filterTree(nodes, nodePaths, pathData);
}