import Application from "./application";
import Executable from "./executable";
import Flavor from "./flavor";
import { type CreateApplicationConfig } from "./tree";
export default class AdeFactory {
    static createApplication({ name, version, build }: CreateApplicationConfig): Application;
    static getExecutables(application: Application): Executable[];
    static getExecutableByName(application: Application, name?: string): Executable;
    static getExecutableByConfig(application: Application, config?: {
        name: string;
    }): Executable;
    static getFlavorsByApplicationVersion(executable: Executable, version: string): Flavor[];
    static getExecutableFlavors(executable: Executable): Flavor[];
    static getFlavorByName(executable: Executable, name?: string): Flavor | undefined;
    static getFlavorByConfig(executable: Executable, config?: {
        name: string;
    }): Flavor | undefined;
}
